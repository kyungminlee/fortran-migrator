#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "symbol_database.h"
#include "source_scanner.h"
#include "transformation.h"
#include "rewriter.h"

static void print_usage(const char* prog) {
    std::cerr
        << "Usage: " << prog << " [options] <source-dir> <output-dir>\n"
        << "\n"
        << "Options:\n"
        << "  --kind <10|16>         Target floating-point KIND (default: 16)\n"
        << "  --library <path>       Static/shared library for symbol extraction\n"
        << "  --dry-run              Show planned transformations without writing\n"
        << "  --convergence-check    Migrate from both S/D (or C/Z) and diff\n"
        << "  -h, --help             Show this help\n";
}

int main(int argc, char* argv[]) {
    int target_kind = 16;
    std::string library_path;
    std::string source_dir;
    std::string output_dir;
    bool dry_run = false;
    bool convergence_check = false;

    // Parse arguments
    std::vector<std::string> positional;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            return 0;
        } else if (arg == "--kind" && i + 1 < argc) {
            target_kind = std::stoi(argv[++i]);
        } else if (arg == "--library" && i + 1 < argc) {
            library_path = argv[++i];
        } else if (arg == "--dry-run") {
            dry_run = true;
        } else if (arg == "--convergence-check") {
            convergence_check = true;
        } else if (arg[0] == '-') {
            std::cerr << "Unknown option: " << arg << "\n";
            return 1;
        } else {
            positional.push_back(arg);
        }
    }

    if (positional.size() < 2) {
        print_usage(argv[0]);
        return 1;
    }

    source_dir = positional[0];
    output_dir = positional[1];

    if (target_kind != 10 && target_kind != 16) {
        std::cerr << "Error: --kind must be 10 or 16\n";
        return 1;
    }

    // TODO: implement pipeline
    std::cout << "fortran-migrator\n"
              << "  source:    " << source_dir << "\n"
              << "  output:    " << output_dir << "\n"
              << "  kind:      " << target_kind << "\n"
              << "  library:   " << (library_path.empty() ? "(none)" : library_path) << "\n"
              << "  dry-run:   " << (dry_run ? "yes" : "no") << "\n"
              << "  convergence-check: " << (convergence_check ? "yes" : "no") << "\n";

    return 0;
}
