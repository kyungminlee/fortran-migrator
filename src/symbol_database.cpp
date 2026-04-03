#include "symbol_database.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <fstream>
#include <regex>
#include <sstream>

namespace fmig {

// Precision prefixes used by BLAS/LAPACK/ScaLAPACK.
// A symbol like "PDGESV" has a 'P' routing prefix and 'D' precision prefix.
static bool is_precision_prefix(char c) {
    return c == 'S' || c == 'D' || c == 'C' || c == 'Z';
}

static std::string to_upper(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return std::toupper(c); });
    return s;
}

char SymbolDatabase::target_prefix(PrecisionKind prec, int kind) {
    switch (prec) {
    case PrecisionKind::SingleReal:
    case PrecisionKind::DoubleReal:
        return (kind == 10) ? 'E' : 'Q';
    case PrecisionKind::SingleComplex:
    case PrecisionKind::DoubleComplex:
        return (kind == 10) ? 'Y' : 'X';
    case PrecisionKind::None:
        return '\0';
    }
    return '\0';
}

void SymbolDatabase::classify_and_insert(const std::string& name) {
    std::string upper = to_upper(name);
    if (upper.empty() || index_.count(upper))
        return;

    SymbolEntry entry;
    entry.original_name = upper;

    // Try to classify the precision prefix.
    // ScaLAPACK routines start with 'P' followed by a precision prefix.
    size_t prefix_pos = 0;
    if (upper.size() > 2 && upper[0] == 'P' && is_precision_prefix(upper[1])) {
        prefix_pos = 1;
    } else if (upper.size() > 1 && is_precision_prefix(upper[0])) {
        prefix_pos = 0;
    } else {
        // No recognized precision prefix.
        entry.base_name = upper;
        entry.precision = PrecisionKind::None;
        index_[upper] = entries_.size();
        entries_.push_back(std::move(entry));
        return;
    }

    char p = upper[prefix_pos];
    switch (p) {
    case 'S': entry.precision = PrecisionKind::SingleReal; break;
    case 'D': entry.precision = PrecisionKind::DoubleReal; break;
    case 'C': entry.precision = PrecisionKind::SingleComplex; break;
    case 'Z': entry.precision = PrecisionKind::DoubleComplex; break;
    default:  entry.precision = PrecisionKind::None; break;
    }

    // base_name is everything except the precision character.
    entry.base_name = upper.substr(0, prefix_pos) + upper.substr(prefix_pos + 1);

    index_[upper] = entries_.size();
    entries_.push_back(std::move(entry));
}

bool SymbolDatabase::load_from_library(const std::filesystem::path& lib_path) {
    // Run nm to extract defined symbols.
    std::string cmd = "nm --defined-only -g \"" + lib_path.string() +
                      "\" 2>/dev/null";
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe)
        return false;

    char buf[1024];
    // nm output: "address type name"
    std::regex nm_re(R"(\s*[0-9a-fA-F]*\s+[TtWw]\s+(\S+))");
    while (fgets(buf, sizeof(buf), pipe)) {
        std::string line(buf);
        std::smatch m;
        if (std::regex_search(line, m, nm_re)) {
            std::string sym = m[1].str();
            // Strip leading underscore (macOS) and trailing underscore (Fortran).
            if (!sym.empty() && sym.front() == '_')
                sym.erase(sym.begin());
            if (!sym.empty() && sym.back() == '_')
                sym.pop_back();
            classify_and_insert(sym);
        }
    }
    pclose(pipe);
    return !entries_.empty();
}

bool SymbolDatabase::load_from_source_dir(const std::filesystem::path& src_dir) {
    // Scan Fortran source files for SUBROUTINE and FUNCTION definitions.
    std::regex def_re(
        R"(^\s{0,5}\s+(SUBROUTINE|FUNCTION|(?:DOUBLE\s+PRECISION\s+FUNCTION)|(?:COMPLEX\*\d+\s+FUNCTION)|(?:REAL\s+FUNCTION)|(?:INTEGER\s+FUNCTION)|(?:LOGICAL\s+FUNCTION))\s+(\w+))",
        std::regex::icase);
    // Simpler pattern for the common case:
    std::regex simple_re(
        R"((?:SUBROUTINE|FUNCTION)\s+([A-Za-z]\w*))",
        std::regex::icase);

    for (auto& entry : std::filesystem::recursive_directory_iterator(src_dir)) {
        if (!entry.is_regular_file())
            continue;
        auto ext = entry.path().extension().string();
        std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
        if (ext != ".f" && ext != ".f90" && ext != ".for" && ext != ".f95")
            continue;

        std::ifstream in(entry.path());
        std::string line;
        while (std::getline(in, line)) {
            // Skip comments (fixed-form: C/c/* in column 1).
            if (!line.empty() &&
                (line[0] == 'C' || line[0] == 'c' || line[0] == '*'))
                continue;

            std::smatch m;
            if (std::regex_search(line, m, simple_re)) {
                classify_and_insert(m[1].str());
            }
        }
    }
    return !entries_.empty();
}

std::string SymbolDatabase::target_name(const std::string& name,
                                        int target_kind) const {
    std::string upper = to_upper(name);
    auto it = index_.find(upper);
    if (it == index_.end())
        return {};

    const auto& entry = entries_[it->second];
    if (entry.precision == PrecisionKind::None)
        return {}; // Type-independent — no rename needed.

    char new_prefix = target_prefix(entry.precision, target_kind);
    if (new_prefix == '\0')
        return {};

    // Reconstruct: insert new prefix at the position where the old one was.
    std::string result = entry.base_name;
    size_t insert_pos = (upper.size() > 2 && upper[0] == 'P' &&
                         is_precision_prefix(upper[1]))
                            ? 1
                            : 0;
    result.insert(result.begin() + insert_pos, new_prefix);
    return result;
}

bool SymbolDatabase::contains(const std::string& name) const {
    return index_.count(to_upper(name)) > 0;
}

} // namespace fmig
