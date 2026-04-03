#pragma once

#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "fortran_format.h"
#include "symbol_database.h"
#include "transformation.h"

namespace fmig {

// Scans a single Fortran source file and produces a list of transformations
// needed to migrate it to the target precision KIND.
class SourceScanner {
public:
    explicit SourceScanner(const SymbolDatabase& symbols, int target_kind);

    // Scan a source file and return all transformations.
    std::vector<Transformation> scan(const std::filesystem::path& file_path) const;

    // Scan source text directly (for testing).
    std::vector<Transformation> scan_text(std::string_view source,
                                          FortranForm form) const;

private:
    // Individual scan passes.
    void scan_type_declarations(std::string_view source, FortranForm form,
                                std::vector<Transformation>& out) const;
    void scan_routine_names(std::string_view source, FortranForm form,
                            std::vector<Transformation>& out) const;
    void scan_call_sites(std::string_view source, FortranForm form,
                         std::vector<Transformation>& out) const;
    void scan_external_decls(std::string_view source, FortranForm form,
                             std::vector<Transformation>& out) const;
    void scan_literals(std::string_view source, FortranForm form,
                       std::vector<Transformation>& out) const;
    void scan_intrinsics(std::string_view source, FortranForm form,
                         std::vector<Transformation>& out) const;
    void scan_xerbla_strings(std::string_view source, FortranForm form,
                             std::vector<Transformation>& out) const;

    const SymbolDatabase& symbols_;
    int target_kind_;
};

} // namespace fmig
