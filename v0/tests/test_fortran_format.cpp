#include "fortran_format.h"

#include <cassert>
#include <iostream>

using namespace fmig;

static void test_detect_form() {
    assert(detect_form(".f") == FortranForm::Fixed);
    assert(detect_form(".for") == FortranForm::Fixed);
    assert(detect_form(".f90") == FortranForm::Free);
    assert(detect_form(".f95") == FortranForm::Free);
    assert(detect_form(".F") == FortranForm::Fixed);
}

static void test_fixed_form_comment() {
    assert(is_fixed_form_comment("C     This is a comment"));
    assert(is_fixed_form_comment("c     this too"));
    assert(is_fixed_form_comment("*     and this"));
    assert(is_fixed_form_comment(""));
    assert(!is_fixed_form_comment("      DOUBLE PRECISION A"));
}

static void test_fixed_form_continuation() {
    assert(!is_fixed_form_continuation("      X = 1"));
    assert(is_fixed_form_continuation("     +    BETA*C"));
    assert(is_fixed_form_continuation("     $    BETA*C"));
    assert(!is_fixed_form_continuation("C     comment"));
}

static void test_fixed_form_statement() {
    std::string_view line = "      DOUBLE PRECISION A(LDA,*),B(LDB,*)";
    auto stmt = fixed_form_statement(line);
    assert(stmt == "DOUBLE PRECISION A(LDA,*),B(LDB,*)");
}

int main() {
    test_detect_form();
    test_fixed_form_comment();
    test_fixed_form_continuation();
    test_fixed_form_statement();
    std::cout << "fortran_format tests passed\n";
    return 0;
}
