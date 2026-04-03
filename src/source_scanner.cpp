#include "source_scanner.h"

#include <fstream>
#include <sstream>

namespace fmig {

SourceScanner::SourceScanner(const SymbolDatabase& symbols, int target_kind)
    : symbols_(symbols), target_kind_(target_kind) {}

std::vector<Transformation>
SourceScanner::scan(const std::filesystem::path& file_path) const {
    std::ifstream in(file_path);
    if (!in)
        return {};

    std::ostringstream buf;
    buf << in.rdbuf();
    std::string source = buf.str();

    FortranForm form = detect_form(file_path.extension().string());
    return scan_text(source, form);
}

std::vector<Transformation>
SourceScanner::scan_text(std::string_view source, FortranForm form) const {
    std::vector<Transformation> transforms;

    scan_type_declarations(source, form, transforms);
    scan_routine_names(source, form, transforms);
    scan_call_sites(source, form, transforms);
    scan_external_decls(source, form, transforms);
    scan_literals(source, form, transforms);
    scan_intrinsics(source, form, transforms);
    scan_xerbla_strings(source, form, transforms);

    validate_and_sort(transforms);
    return transforms;
}

// --- Stub implementations (to be filled in) ---

void SourceScanner::scan_type_declarations(
    std::string_view /*source*/, FortranForm /*form*/,
    std::vector<Transformation>& /*out*/) const {
    // TODO: Identify REAL, DOUBLE PRECISION, COMPLEX, DOUBLE COMPLEX,
    //       REAL*8, COMPLEX*16, etc. and replace with REAL(KIND=k) or
    //       COMPLEX(KIND=k).
}

void SourceScanner::scan_routine_names(
    std::string_view /*source*/, FortranForm /*form*/,
    std::vector<Transformation>& /*out*/) const {
    // TODO: Find SUBROUTINE/FUNCTION/ENTRY definition lines and rename
    //       using symbol database.
}

void SourceScanner::scan_call_sites(
    std::string_view /*source*/, FortranForm /*form*/,
    std::vector<Transformation>& /*out*/) const {
    // TODO: Find CALL statements and function references in expressions,
    //       rename using symbol database.
}

void SourceScanner::scan_external_decls(
    std::string_view /*source*/, FortranForm /*form*/,
    std::vector<Transformation>& /*out*/) const {
    // TODO: Find EXTERNAL declarations and rename listed symbols.
}

void SourceScanner::scan_literals(
    std::string_view /*source*/, FortranForm /*form*/,
    std::vector<Transformation>& /*out*/) const {
    // TODO: Find floating-point literals with D/E exponents and convert
    //       to target KIND suffix form.
}

void SourceScanner::scan_intrinsics(
    std::string_view /*source*/, FortranForm /*form*/,
    std::vector<Transformation>& /*out*/) const {
    // TODO: Find type-specific intrinsic calls (DBLE, DCMPLX, DCONJG, etc.)
    //       and convert to generic equivalents with KIND argument.
}

void SourceScanner::scan_xerbla_strings(
    std::string_view /*source*/, FortranForm /*form*/,
    std::vector<Transformation>& /*out*/) const {
    // TODO: Find string literals in XERBLA calls containing routine names
    //       and update the prefix.
}

} // namespace fmig
