#include "rewriter.h"

#include <cassert>
#include <iostream>

using namespace fmig;

static void test_target_filename() {
    // Without a populated symbol database, target_filename returns the
    // original name.  Full rename tests require a loaded database.
    SymbolDatabase db;
    std::string result = target_filename("dgemm.f", db, 16);
    // DB is empty, so no rename.
    assert(result == "dgemm.f");
    std::cout << "test_target_filename: basic case passed\n";
}

int main() {
    test_target_filename();
    std::cout << "rewriter tests passed\n";
    return 0;
}
