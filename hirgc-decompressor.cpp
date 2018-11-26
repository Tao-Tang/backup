# include <cstdlib>
# include <iostream>
# include <fstream>
# include <ctime>
# include <sys/time.h>
# include <cmath>
# include <vector>
# include <algorithm>
# include <map>
# include <climits>
# include <cstdint>
# include <cstring>
# include <string>
# include <utility>
# include <unistd.h>

using namespace std;

# include "util.h"

//for target and ref genomic
char *meta_data;
POSITION_RANGE *pos_vec, *n_vec;
POSITION_OTHER_CHAR *other_char_vec;
int *line_break_vec;
char *ref_seq;
char *tar_seq;
int other_char_len, str_len, tar_seq_len, ref_seq_len, pos_vec_len, n_vec_len, line_break_len;



void initial() 
{
	meta_data = new char[2048];
	tar_seq = new char[MAX_CHAR_NUM];
	ref_seq = new char[MAX_CHAR_NUM];
	return;
}

void clear() {
	delete[] meta_data;
	delete[] tar_seq;
	delete[] ref_seq;
}

void readRefFile(char *refFile) {
	ref_seq_len = 0;
	char ch[2048];
	FILE *fp = fopen(refFile, "r");
	if (NULL == fp) {
		printf("fail to open file %s\n", refFile);
		exit(0);
		return;
	}
	int temp_len;
	char temp_ch;
	fgets(ch, 2048, fp);

	while (fscanf(fp, "%s", ch) != EOF) {
		// printf("%s\n", ch);
		temp_len = strlen(ch);
		for (int i = 0; i < temp_len; i++) {
			temp_ch = ch[i];
			if (islower(temp_ch)) {
				temp_ch = toupper(temp_ch);
			}
			if (temp_ch == 'A' || temp_ch == 'G' || temp_ch == 'C' || temp_ch == 'T' ) {
				ref_seq[ref_seq_len++] = temp_ch;
			}
		}
	}
	ref_seq[ref_seq_len] = '\0';
	fclose(fp);
}

void readRunLengthCoding(FILE *fp, int &vec_len, int **vec) {
	int code_len, temp;
	vector<int> code;
	fscanf(fp, "%d", &code_len);
	for (int i = 0; i < code_len; i++) {
		fscanf(fp, "%d", &temp);
		code.push_back(temp);
	}

	vec_len = 0;
	for (int i = 1; i < code_len; i += 2) {
		vec_len += code[i];
	}

	if (vec_len > 0) {
		*vec = new int[vec_len];
		int k = 0;
		for (int i = 0; i < code_len; i += 2) {
			for (int j = 0; j < code[i+1]; j++) {
				(*vec)[k++] = code[i];
			}
		}
	}
}

void readOtherData(FILE *fp) {
	meta_data = new char[2048];
	fgets(meta_data, 2048, fp);
	// fscanf(fp, "%s", meta_data);

	readRunLengthCoding(fp, line_break_len, &line_break_vec);

	fscanf(fp, "%d", &pos_vec_len);
	// printf("pos_vec_len: %d\n", pos_vec_len);
	if (pos_vec_len > 0) {
		pos_vec = new POSITION_RANGE[pos_vec_len];
		for (int i = 0; i < pos_vec_len; i++) {
			fscanf(fp, "%d%d", &pos_vec[i].begin, &pos_vec[i].length);
		}
	}

	fscanf(fp, "%d", &n_vec_len);
	// printf("n_vec_len: %d\n", n_vec_len);
	if (n_vec_len > 0) {
		n_vec = new POSITION_RANGE[n_vec_len];
		for (int i = 0; i < n_vec_len; i++) {
			fscanf(fp, "%d%d", &n_vec[i].begin, &n_vec[i].length);
		}
	}

	fscanf(fp, "%d", &other_char_len);
	// printf("other_char_len: %d\n", other_char_len);
	if (other_char_len > 0) {
		other_char_vec = new POSITION_OTHER_CHAR[other_char_len];
		for (int i = 0; i < other_char_len; i++) {
			fscanf(fp, "%d", &other_char_vec[i].pos);
		}
		//read other_char_vec information
		vector<int> arr;
		int size, temp;
		fscanf(fp, "%d", &size);
		for (int i = 0; i < size; i++) {
			fscanf(fp, "%d", &temp);
			arr.push_back(temp);
		}
		if (size < 10) {
			char *temp_str = new char[other_char_len+10];
			fscanf(fp, "%s", temp_str);
			for (int i = 0; i < other_char_len; i++) {
				other_char_vec[i].ch = arr[temp_str[i]-'0'];
			}
			delete[] temp_str;
		} else {
			int bit_num = ceil(log(size)/log(2));
			int v_num = floor(32.0/bit_num);

			for (int i = 0; i < other_char_len; ) {
				unsigned int v;
				fscanf(fp, "%u", &v);
				vector<int> temp_arr;

				int temp_i = i;
				for (int j = 0; j < v_num && temp_i < other_char_len; j++, temp_i++) {
					int mod = v%(1<<bit_num);
					v >>= bit_num;
					temp_arr.push_back(arr[mod]);
				}
				for (int j = temp_arr.size()-1; j >= 0 && i < other_char_len; j--, i++) {
					other_char_vec[i].ch = temp_arr[j];
				}
			}
		}
	}
}

bool exitSpace(char *str) {
	for (int i = 0; str[i]; i++) {
		if (str[i] == ' ') {
			return true;
		}
	}
	return false;
} 

void readCompressedFile(char *tarFile) {
	FILE *fp = fopen(tarFile, "r");
	if (NULL == fp) {
		printf("ERROR! fail to open file %s\n", tarFile);
		exit(0);
	}
	readOtherData(fp);
	// cerr << "after readOtherData()...\n";
	int pre_pos = 0, cur_pos, length;;
	tar_seq_len = 0;
	// printf("%s\n", ref_seq);
	char *str = new char[MAX_CHAR_NUM];

	while (fgets(str, MAX_CHAR_NUM, fp) != NULL) {
		if (str[0] == '\n') continue;
		// printf("cmd: %s", str);
		// cerr << "str: " << str << "\n";
		if (exitSpace(str)) {
			sscanf(str, "%d%d", &cur_pos, &length);
			// cerr << cur_pos << "; " << length << "\n";
			cur_pos += pre_pos;
			length += sub_str_num;
			pre_pos = cur_pos + length;

			// cerr << cur_pos << "; " << length << "\n";

			for (int i = cur_pos, j = 0; j < length; j++, i++) {
				tar_seq[tar_seq_len++] = ref_seq[i];
			}
			tar_seq[tar_seq_len] = '\0';
			// cerr<<"over\n";
		} else {
			int str_len = strlen(str);
			// cerr << str_len << "\n";
			for (int i = 0; i < str_len-1; i++) {
				tar_seq[tar_seq_len++] = invert_code_rule[str[i]-'0'];
			}
		}
	}
	tar_seq[tar_seq_len] = '\0';
	// printf("%s\n", tar_seq);
	fclose(fp);
	delete[] str;
}
 
void saveDecompressedData(char *resultFile) {
	FILE *fp = fopen(resultFile, "w");
	if (NULL == fp) {
		printf("ERROR! fail to open file %s\n", resultFile);
		exit(0);
		return;
	}

	for (int i = 1; i < other_char_len; i++) {
		other_char_vec[i].pos += other_char_vec[i-1].pos;
	}

	char *temp_seq = new char[MAX_CHAR_NUM];
	
	strcpy(temp_seq, tar_seq);
	int tt = 0, j = 0;

	for (int i = 0; i < other_char_len; i++) {
		while (tt < other_char_vec[i].pos && tt < tar_seq_len) {
			tar_seq[j++] = temp_seq[tt++];
		}
		tar_seq[j++] = other_char_vec[i].ch + 'A';
	}
	while (tt < tar_seq_len) {
		tar_seq[j++] = temp_seq[tt++];
	}
	tar_seq[j] = '\0';
	tar_seq_len = j;

	int str_len = 0;
	int r = 0; 

	char *str = new char[MAX_CHAR_NUM];

	for (int i = 0; i < n_vec_len; i++) {
		for (int j = 0; j < n_vec[i].begin; j++) {
			str[str_len++] = tar_seq[r++];
		}
		for (int j = 0; j < n_vec[i].length; j++) {
			str[str_len++] = 'N';
		}
	}
	while (r < tar_seq_len) {
		str[str_len++] = tar_seq[r++];
	}
	str[str_len] = '\0';

	fprintf(fp, "%s", meta_data);

	int k = 0;
	for (int i = 0; i < pos_vec_len; i++) {
		k += pos_vec[i].begin;
		int temp = pos_vec[i].length;
		for (int j = 0; j < temp; j++) {
			str[k] = tolower(str[k]);
			k++;
		}
	}
	
	for (int i = 1; i < line_break_len; i++) {
		line_break_vec[i] += line_break_vec[i-1];
	}
	int k_lb = 0;
	
	int temp_seq_len = 0;
	for (int i = 0; i < str_len; i++) {
		if(i == line_break_vec[k_lb]) {
			temp_seq[temp_seq_len++] = '\n';
			k_lb++;
		}
		temp_seq[temp_seq_len++] = str[i];
	}
	temp_seq[temp_seq_len] = '\0';
	fprintf(fp, "%s", temp_seq);
	// printf("k_lb: %d\n", k_lb);
	while ((k_lb++) < line_break_len) {
		fprintf(fp, "\n");
	}
	fclose(fp);

	delete[] temp_seq;
	delete[] str;

	if (pos_vec_len > 0) delete[] pos_vec;
	delete[] line_break_vec;
	if (n_vec_len > 0) delete[] n_vec;
	if (other_char_len > 0) delete[] other_char_vec;
}

void decompressFile(char *refFile, char *tarFile, char *resultFile) {
	readRefFile(refFile);//must read first;
	// cerr << "after readRefFile()...\n";
	readCompressedFile(tarFile);
	// cerr << "after readCompressedFile()...\n";
	saveDecompressedData(resultFile);
}

void decompressFile(char *refFile, char *tarFile) {
	initial();//****important
	char cmd[1024]; 
	sprintf(cmd, "7za x %s", tarFile);
	system(cmd);

	printf("decompressing...\n");
	
	char temp[1024];
	sprintf(temp, "%s", tarFile);

	char res[1024];
	sprintf(res, "dec_%s", temp);
	decompressFile(refFile, temp, res);

	clear();
}


bool getList(char *list_file, vector<string> &name_list) {
	FILE *fp = fopen(list_file, "r");
	if (fp == NULL) {
		printf("%s open fail!\n", list_file);
		return false;
	}
	char str[100];
	while (fscanf(fp, "%s", str) != EOF) {
		name_list.push_back(string(str));
	}
	fclose(fp);
	if (name_list.size() == 0) {
		printf("%s is empty!\n", list_file);
		return false;
	}
	return true;
}

class deCompressor
{
public:
	string ref_name;
	deCompressor(string m_ref_name)
	{
		ref_name = m_ref_name;
		initial();
		char buffer[1024];
		sprintf(buffer,"%s",ref_name.c_str());
		readRefFile(buffer);
	}

	void decompress(char *new_ref, char *target)
	{
		if ( ref_name != string(new_ref)  )
		{
			readRefFile(new_ref);
			ref_name = string(new_ref);
		}

		char src[1024];
		sprintf(src,"%s_ref_%s",target,new_ref);
		printf("%s\n",target);
		readCompressedFile(src);
	// cerr << "after readCompressedFile()...\n";
		saveDecompressedData(target);
		char cmd[1024];
		sprintf(cmd,"rm %s",src);
		system(cmd);
		return;
	}
	~deCompressor() { clear(); }

};
