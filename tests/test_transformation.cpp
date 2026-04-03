#include "transformation.h"

#include <cassert>
#include <iostream>

using namespace fmig;

static void test_no_overlap() {
    std::vector<Transformation> ts = {
        {10, 5, "HELLO", "first"},
        {20, 3, "BYE", "second"},
    };
    assert(validate_and_sort(ts));
    // Should be sorted descending by offset.
    assert(ts[0].offset == 20);
    assert(ts[1].offset == 10);
}

static void test_overlap_detected() {
    std::vector<Transformation> ts = {
        {10, 10, "AAA", "first"},  // covers [10, 20)
        {15, 5, "BBB", "second"},  // covers [15, 20) — overlaps
    };
    assert(!validate_and_sort(ts));
}

static void test_apply() {
    std::string source = "DOUBLE PRECISION A, B";
    std::vector<Transformation> ts = {
        {0, 16, "REAL(KIND=16)", "type"},
    };
    validate_and_sort(ts);
    std::string result = apply_transformations(source, ts);
    assert(result == "REAL(KIND=16) A, B");
}

static void test_apply_multiple() {
    std::string source = "CALL DGEMM(...) ! uses DGEMM";
    std::vector<Transformation> ts = {
        {5, 5, "QGEMM", "call rename"},
        {23, 5, "QGEMM", "comment rename"},
    };
    validate_and_sort(ts);
    std::string result = apply_transformations(source, ts);
    assert(result == "CALL QGEMM(...) ! uses QGEMM");
}

int main() {
    test_no_overlap();
    test_overlap_detected();
    test_apply();
    test_apply_multiple();
    std::cout << "transformation tests passed\n";
    return 0;
}
