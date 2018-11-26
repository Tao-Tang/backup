#ifndef UTIL_H
#define UTIL_H
const int MAX_CHAR_NUM = 1<<28;
// const unsigned char code_rule[4] = {0, 1, 2, 3};//A-0; C-1; G-2; T-3;
const char invert_code_rule[4] = {'A', 'C', 'G', 'T'};
const int sub_str_num = 20;
struct POSITION_RANGE { //record lower case range
	int begin, length;//
};

struct POSITION_OTHER_CHAR {
	int pos, ch;
};

struct REF_NODE {
	int v, next;
};
#endif