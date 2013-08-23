#pragma once
#ifndef TEST_H__INCLUDED
#define TEST_H__INCLUDED

// assertions for "unit testing"
void assertEqual(int expected, int actual, const char* message);
void assertEqual_u8(uint8_t expected, uint8_t actual, const char* message);
bool is_odd(int x);
bool is_even(int x);

#endif