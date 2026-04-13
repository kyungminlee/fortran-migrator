#pragma once

#include <filesystem>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace fmig {

// Classification of a routine's precision prefix.
enum class PrecisionKind {
    SingleReal,    // S prefix
    DoubleReal,    // D prefix
    SingleComplex, // C prefix
    DoubleComplex, // Z prefix
    None,          // No precision prefix (e.g., LSAME, XERBLA)
};

// A routine entry in the symbol database.
struct SymbolEntry {
    std::string original_name; // e.g., "DGEMM"
    std::string base_name;     // e.g., "GEMM"
    PrecisionKind precision;   // e.g., DoubleReal
};

// Maps known library symbols to their precision classification and provides
// rename lookups for a given target KIND.
class SymbolDatabase {
public:
    // Build from a compiled library archive (.a) or shared library (.so/.dylib)
    // by running `nm` and parsing the output.
    bool load_from_library(const std::filesystem::path& lib_path);

    // Build from a source directory by scanning for SUBROUTINE/FUNCTION
    // definitions in Fortran source files.
    bool load_from_source_dir(const std::filesystem::path& src_dir);

    // Look up the target name for a given routine at a specific target KIND.
    // Returns empty string if the symbol is not in the database or has no
    // precision prefix (type-independent routine).
    std::string target_name(const std::string& name, int target_kind) const;

    // Check whether a symbol is a known library routine.
    bool contains(const std::string& name) const;

    // All entries.
    const std::vector<SymbolEntry>& entries() const { return entries_; }

private:
    void classify_and_insert(const std::string& name);

    // Returns the new prefix character for a given precision and target KIND.
    static char target_prefix(PrecisionKind prec, int kind);

    std::vector<SymbolEntry> entries_;
    std::unordered_map<std::string, size_t> index_; // name → entries_ index
};

} // namespace fmig
