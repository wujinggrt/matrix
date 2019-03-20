#ifndef UNIT_TEST_H__
#define UNIT_TEST_H__

#include <iostream>
#include <string>

static int test_count = 0;
static int test_pass = 0;

#define TEST_CASE(NAME)  \
    extern void TEST_CASE_FUNCTION_##NAME(); \
    class TEST_CASE_CLASS_##NAME \
    { \
    public: \
        TEST_CASE_CLASS_##NAME() \
        { \
            std::cout << #NAME << "\n"; \
            std::cout << __TIME__ << "\n"; \
            std::cout << __LINE__ << " " << __FILE__ << "\n"; \
            TEST_CASE_FUNCTION_##NAME(); \
            std::cout << "\n"; \
        } \
    } TEST_CASE_INSTANCE_##NAME; \
    void TEST_CASE_FUNCTION_##NAME()

#define EXPECT_TRUE(condition) do{ if (!(condition)) throw 0; } while(0)

#define TEST_PRINT(MESSAGE) std::cout << (MESSAGE) << std::endl;

#endif