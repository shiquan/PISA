/*  
    Copyright (C) 2016,2017  BGI Research

    Author: Shi Quan (shiquan@genomics.cn)

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
    THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE. 
*/

#ifndef NUMBER_HEADER
#define NUMBER_HEADER
#include <stdio.h>
#include <stdlib.h>

extern int get_numbase(const char *s);
extern int get_numbase_l(const char *s, int l);
extern int is_ieee_magic_val(const char *val);
extern double nondec2num(char *str, int length);
extern int check_num_likely(const char *str);
extern int check_num_likely_l(const char *str, int length);
extern int check_char_num(const char x);
extern double force2num(char *str);
extern double force2num_l(char *str, int l);
extern int str2int(const char *str);
extern int str2int_l(const char *str, int l);
extern int human2int(const char *str);
#endif
