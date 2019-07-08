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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include <ctype.h>
#include "htslib/kstring.h"
#include <zlib.h>

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_GRAY    "\x1b[37m"
#define ANSI_COLOR_RESET   "\x1b[0m"

#define KSTRING_INIT {0, 0, 0}
enum enclose_type {
    braces = -1,
    bracket,
    rev_braces,
    rev_bracket,
    quote,
    double_quote,
};
struct node {
    int start;
    int pos;
    enum enclose_type type;
    int line;
};
struct stack {
    int m, l;
    struct node *nodes;
};

static struct stack stack = {0, 0, 0};
static int is_paired = 1;
static int is_key = 0;

static int stack_destroy()
{
    free(stack.nodes);
    return 0;
}
static int complement_error(kstring_t *string)
{
    assert (stack.l > 0);
    int i, j = 0;
    i = stack.l == 1 ? 0 : stack.l -2;
        
    for (; i < stack.l; ++i) {
        struct node *node = &stack.nodes[i];
        for (; j < node->start; ++j)
            fprintf(stderr, ANSI_COLOR_GRAY "%c" ANSI_COLOR_RESET, string->s[j]);
        
        for (j = node->start; j < node->pos; ++j) 
            fprintf(stderr, "%c", string->s[j]);
        
        fprintf(stderr, ANSI_COLOR_RED "%c" ANSI_COLOR_RESET, string->s[node->pos]);

        for (j = node->pos + 1; j < string->l && string->s[j] != '\n'; ++j)
            fprintf(stderr, "%c", string->s[j]);
        fprintf(stderr, "\n");
        fprintf(stderr, ANSI_COLOR_GREEN "%*s  line : %d\n" ANSI_COLOR_RESET, node->pos - node->start + 1, "~", node->line);
        j++;
    }
    return 0;
}
static int format_error(kstring_t *string, int start, int pos, int line)
{
    int i;
    for (i = 0; i < start; ++i) {
        fprintf(stderr, ANSI_COLOR_GRAY "%c" ANSI_COLOR_RESET, string->s[i]);
    }
    for (i = start; i < pos; ++i) {
        fprintf(stderr, "%c", string->s[i]);
    }
    for (i = pos; i < string->l && string->s[i] != '\n' && string->s[i] != '{' && string->s[i] != '"' &&
              string->s[i] != '\'' && string->s[i] != '['; ++i)
        fprintf(stderr, ANSI_COLOR_RED "%c" ANSI_COLOR_RESET, string->s[i]);
    fprintf(stderr, "\n");
    fprintf(stderr, ANSI_COLOR_GREEN "%*s line : %d" ANSI_COLOR_RESET "\n", pos - start + 1, "~", line);             
    return 0;
}
static int parse_comment_line(kstring_t *string)
{
    if (string->l == 0)
        return 0;

    char *ss = string->s;
    char *se = string->s + string->l -1;
    while (ss && isspace(*ss))
        ss++;

    if (ss > se || *ss == '#') {
        string->l = 0;
        return 0;
    }

    while (se && isspace(*se))
        se--;

    if (ss != string->s || se - ss + 1 != string->l) {
        string->l = se -ss + 1;
        memmove(string->s, ss, string->l);
        string->s[string->l] = '\0';
    }

    if (string->l == 1)
        return 1;
    
    ss = string->s;
    se = string->s + string->l - 1;
    int mark = 0;

    for (;;) {
        if (ss == se)
            break;
        if (*ss == '/') {
            if (mark == 1) {
                string->l = ss - string->s - 1;
                string->s[string->l] = '\0';
                return string->l ? 1 : 0;    
            } else {
                mark = 1;
            }
        } else {
            mark = 0;
        }
        ss++;
    }
    return string->l ? 1 : 0;
}
static int stack_push(int c, int start, int pos, int line)
{
    if (stack.l == stack.m) {
        stack.m += 4;
        stack.nodes = (struct node *)realloc(stack.nodes, stack.m * sizeof(struct node));
    }    
    assert(start <= pos);
    struct node *node = &stack.nodes[stack.l++];
    node->start = start;
    node->pos = pos;
    node->line = line;
    if (c== '{') {
        node->type = braces;
    } else if (c== '[') {
        node->type = bracket;
    } else if (c == '"') {
        node->type = double_quote;
    } else if (c == '\'') {
        node->type = quote;
    } else if (c == '}') {
        node->type = rev_braces;        
    } else if (c == ']') {
        node->type = rev_bracket;        
    } else {
        fprintf(stderr, "Unknown character, push failed. %c", c);
        return 0;
    }
    if (stack.l > 1) {
        switch (stack.nodes[stack.l-1].type) {
            case double_quote :
                if (stack.nodes[stack.l -2].type == double_quote) {
                    stack.l -= 2;
                    if (is_key == 0) {
                        is_key = 1;
                        is_paired = 0;
                    } else {
                        is_key = 0;
                        is_paired = 1;
                    }                        
                } 
                break;
                
            case quote :
                if (stack.nodes[stack.l - 2].type == quote) {
                    stack.l -= 2;
                    if (is_key == 0) {
                        is_key = 1;
                        is_paired = 0;
                    } else {
                        is_paired = 1;
                        is_key = 0;
                    }
                }
                break;

            case rev_braces :
                if (stack.nodes[stack.l - 2].type == braces) {
                    stack.l -= 2;
                } else {
                    return 1;
                }
                break;

            case rev_bracket :
                if (stack.nodes[stack.l - 2].type == bracket) {
                    stack.l -= 2;
                } else {
                    return 1;
                }
                break;

            case braces :
            case bracket :
                if (stack.nodes[stack.l -2].type == double_quote ||  stack.nodes[stack.l -2].type == quote) {
                    return 1;
                } else {
                    if (is_key == 1) { // here just check the depthest pairs
                        is_key = 0;
                        is_paired = 0;
                    }
                }
                break;

            default :
                fprintf(stderr, "Unknown complenment type.");
                break;
        }
    }
    else {
        if (stack.nodes[0].type == rev_braces || stack.nodes[0].type == rev_bracket)
            return 1;
    }
    
    return 0;
}
static int stack_pop(int c, int start, int pos, int line)
{
    return stack_push(c, start, pos, line);
}
static int stack_quote(int c, int start, int pos, int line)
{
    return stack_push(c, start, pos, line);
}
static int stack_check_quote()
{
    if (stack.l == 0)
        return 0;
    return stack.nodes[stack.l -1].type == double_quote || stack.nodes[stack.l -1].type == quote;
}
static int check_key_pair()
{
    return is_paired;
}
static int check_complement(kstring_t *string)
{
    int i;
    int start = 0;    
    int line = 1;
    int is_val = 0;
    for (i = 0; i < string->l; ++i) {

        switch(string->s[i]) {
            case '\n':
                start = i + 1, line++;
                break;
            case '{' :
            case '[' :
                if (is_val == 0) {
                    if (stack_push(string->s[i], start, i, line)) {
                        complement_error(string);
                        return 1;
                    }
                }
                break;
                
            case '}' :
            case ']' :
                if (is_val == 0) {
                    if (stack_pop(string->s[i], start, i, line)) {
                        complement_error(string);
                        return 1;
                    }
                }
                break;

            case '"' :
            case '\'':
                if (stack_quote(string->s[i], start, i, line)) {
                    complement_error(string);
                    return 1;
                }
                else {
                    is_val = is_val ? 0 : 1;
                }
                
                break;

            case ':' :
                // check the key
                if (stack_check_quote() == 0) { 
                    if (check_key_pair() == 1) {
                        format_error(string, start, i, line);
                        return 1;
                    }
                }
                break;
                
            case ',' :
                // check the key value pair
                /* if (stack_check_quote() == 0) { */
                /*     if (check_key_pair() == 0) { */
                /*         format_error(string, start, i, line); */
                /*         return 1; */
                /*     } */
                /* } */
                /* break; */

            case '\t':
            case ' ': // blank in key or value string
                break;
                
            default :
                if (stack.l == 0 || stack_check_quote() == 0) {
                    format_error(string, start, i, line);
                    return 1;
                }
                break;
        }        
    }
    stack_destroy();
    return 0;
}
char *json_config_open(const char *fname)
{
    FILE *fp;
    fp = fopen(fname, "r");
    if (fp == NULL) {
        fprintf(stderr, "%s : %s\n", fname, strerror(errno));
        return NULL;
    }
    kstring_t string = KSTRING_INIT;
    kstring_t temp = KSTRING_INIT;
    while (kgetline(&temp, fgets, fp) >= 0) {
        if (parse_comment_line(&temp)) {
            kputs(temp.s, &string);
            kputc('\n', &string);
        }
        temp.l = 0;
    }
    fclose(fp);
    if (temp.m)
        free(temp.s);
    //fprintf(stderr, "%s\n", string.s);
    if (check_complement(&string))
        return NULL;
    return string.s;
}

#ifdef JSON_CONFIG_MAIN
int main(int argc, char **argv)
{
    if (argc == 1) {
        fprintf(stderr, "%s   config.json\n", argv[0]);
        return 1;
    }
    char *string = json_config_open(argv[1]);
    if (string) {
        printf("%s", string);
        free(string);
    } else {
        fprintf(stderr, "Failed to parse %s.\n", argv[1]);        
    }
    return 0;
}
#endif
