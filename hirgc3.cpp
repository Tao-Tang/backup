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
# include <thread>
# include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

using namespace std;
 
const int MAX_CHAR_NUM = 1<<28;//maximum length of a chromosome
const int code_rule[4] = {0, 1, 2, 3};//A-0; C-1; G-2; T-3; encoding rule //discarded
const char invert_code_rule[4] = {'A', 'C', 'G', 'T'}; //decoding rule
const int max_arr_num_bit = 30; //relate to length of hash table
const int max_arr_num_bit_shift = max_arr_num_bit>>1; //half
const int max_arr_num = 1<<max_arr_num_bit; // maximum length of hash table
const int min_size = 1<<20; //minimum length for other arrays 

const int sub_str_num = 20; //hash function; the modulus s

struct POSITION_RANGE { //an interval
	int begin, length;//
};

struct POSITION_OTHER_CHAR {//for other characters
	int pos, ch;
};

// this is set to compress the given range of tasks of a file list
class hirgc
{

	public:	
		char *meta_data; // identifier
		POSITION_RANGE *pos_vec, *n_vec; //pos_ver is intervals of low-case letters; n_vec for the letter 'N'
		POSITION_OTHER_CHAR *other_char_vec; // for other characters
		int *line_break_vec;//EOL character; equivalent length of short sequence
		int *tar_seq_code, *ref_seq_code; // target and reference sequence after encoding
		int other_char_len, str_len, tar_seq_len, ref_seq_len, pos_vec_len, n_vec_len, line_break_len; //lengths related to their arrays
		int *loc; //an array of header points
		int *point; // an array of entries
		char *dismatched_str; //mismatched subsequence

	vector<string> genome_list;
	vector<int> cluster_list;
	vector<int> centroid_list;
	int start,end,tid;

	struct  timeval t_start;
	struct  timeval t_end;
	unsigned long timer;	

	hirgc(vector<string> m_genome_list,int m_start, int m_end,int m_tid,vector<int>m_cluster_list={},vector<int> m_centroid_list={}): 
	genome_list(m_genome_list),cluster_list(m_cluster_list), centroid_list(m_centroid_list), start(m_start), end(m_end), tid(m_tid)
	{

	}

	void initial() 
	{ // malloc momories
		meta_data = new char[1024];
		pos_vec = new POSITION_RANGE[min_size];
		line_break_vec = new int[1<<23];
		n_vec = new  POSITION_RANGE[min_size];
		other_char_vec = new POSITION_OTHER_CHAR[min_size];

		tar_seq_code = new int[MAX_CHAR_NUM];
		ref_seq_code = new int[MAX_CHAR_NUM];

		loc = new int[MAX_CHAR_NUM];
		point = new int[max_arr_num];

		dismatched_str = new char[min_size];
	}

void clear() { // free
//	printf("clear 1\n");
	delete[] meta_data;
//	printf("clear 2\n");	
	delete[] pos_vec;
//	printf("clear 3\n");	
	delete[] line_break_vec;
//	printf("clear 4\n");	
	delete[] n_vec;
//	printf("clear 5\n");	
	delete[] other_char_vec;
//	printf("clear 6\n");	
	delete[] tar_seq_code;
//	printf("clear 7\n");	
	delete[] ref_seq_code;
//	printf("clear 8\n");	
	delete[] loc;
//	printf("clear 9\n");	
	delete[] point;
	printf("clear finish\n");	
//	delete[] dismatched_str;
}

	virtual void readRefFile(char *refFile);
	virtual void readTarFile(char *tarFile);
	virtual void writeRunLengthCoding(FILE *fp, int vec_len, int *vec);
	virtual void saveOtherData(FILE *fp);
	virtual void preProcessRef();
	virtual void searchMatch(char *refFile);
	virtual void hirgc_cluster_compress(char *res_folder);
	virtual void hirgc_single_ref_compress(char *res_folder, char *ref);
	virtual void write_time_record(char *res_folder);


};



int acgtTrans(char ch) 
{ //encoding rule
    if (ch == 'A') {
        return 0;
    }
    if (ch == 'C') {
        return 1;
    } 
    if (ch == 'G') {
        return 2;
    }
    if (ch == 'T') {
        return 3;
    }
    return 4;
}

string strip_string(string s)
{
	string result;
	if(s.find("/") != string::npos)
	{
		result = string( s.substr(s.find_last_of("/")+1 ) );
		return result;
	}
	return s;
}



void hirgc::readRefFile(char *refFile) 
{ // processing reference file
	int _ref_seq_len = 0;
	char ch[1024];
	FILE *fp = fopen(refFile, "r");
	if (NULL == fp) {
		printf("fail to open file %s\n", refFile);
		return;
	}
	int temp_len, index;
	char temp_ch;
	
	// fscanf(fp, "%s", ch);//meta_data
	fgets(ch, 1024, fp);
	while (fscanf(fp, "%s", ch) != EOF) {
		temp_len = strlen(ch);
		for (int i = 0; i < temp_len; i++) {
			temp_ch = ch[i];
			if (islower(temp_ch)) {
				temp_ch = toupper(temp_ch);
			}
			index = acgtTrans(temp_ch);
			if (index^4) {//only A T C T are saved
				ref_seq_code[_ref_seq_len++] = index;
			}
		}
	}
	fclose(fp);
	ref_seq_len = _ref_seq_len;
}


void hirgc::readTarFile(char *tarFile) {// processing target file; recording all auxiliary information  
	
	FILE *fp = fopen(tarFile, "r");
	if (NULL == fp) {
		printf("fail to open file %s\n", tarFile);
		return;
	}
	char ch[1024], chr;
	
	int _tar_seq_len = 0;
	other_char_len = 0, n_vec_len = 0;
	pos_vec_len = line_break_len = 0;

	int letters_len = 0, n_letters_len = 0, index, ch_len;
	bool flag = true, n_flag = false; // first is upper case //first is not N

	// fscanf(fp, "%s", meta_data); //meta_data
	fgets(meta_data, 1024, fp);

	while (fscanf(fp, "%s", ch) != EOF) {
		ch_len = strlen(ch);

		for (int i = 0; i < ch_len; i++) {
			chr = ch[i];
			if (islower(chr)) {
				if (flag) { //previous is upper case
					flag = false; //change status of switch
					pos_vec[pos_vec_len].begin = (letters_len);
					letters_len = 0;
				}
				chr = toupper(chr);
			} else {
				if (isupper(chr)) {
					if (!flag) {
						flag = true;
						pos_vec[pos_vec_len].length = (letters_len);
						pos_vec_len++;
						letters_len = 0;
					}
				}
			}
			letters_len++;

			//ch is an upper letter
			if (chr != 'N') {
				index = acgtTrans(chr);
				if (index^4) {
					// tar_seq_code[tar_seq_len++] = code_rule[index];
					tar_seq_code[_tar_seq_len++] = index;
				} else {
					other_char_vec[other_char_len].pos = _tar_seq_len;
					other_char_vec[other_char_len].ch = chr-'A';
					other_char_len++;
				}
			}

			if (!n_flag) {
				if (chr == 'N') {
					n_vec[n_vec_len].begin = n_letters_len;
					n_letters_len = 0;
					n_flag = true;
				}
			} else {//n_flag = true
				if (chr != 'N'){
					n_vec[n_vec_len].length = n_letters_len;
					n_vec_len++;
					n_letters_len = 0;
					n_flag = false;
				}
			}
			n_letters_len++;

		}
		line_break_vec[line_break_len++] = ch_len;
	}

	if (!flag) {
		pos_vec[pos_vec_len].length = (letters_len);
		pos_vec_len++;
	}
	
	if (n_flag) {
		n_vec[n_vec_len].length = (n_letters_len);
		n_vec_len++;
	}

	for (int i = other_char_len-1; i > 0; i--) {
		other_char_vec[i].pos -= other_char_vec[i-1].pos;
	}
	fclose(fp);
	tar_seq_len = _tar_seq_len;
}

void hirgc::writeRunLengthCoding(FILE *fp, int vec_len, int *vec) { // run-length coding for EOL
	vector<int> code;
	if (vec_len > 0) {
		code.push_back(vec[0]);
		int pre_value = code[0], cnt = 1;
		for (int i = 1; i < vec_len; i++) {
			if (vec[i] == pre_value) {
				cnt++;
			} else {
				code.push_back(cnt);
				code.push_back(vec[i]);
				pre_value = vec[i];
				cnt = 1;
			}
		}
		code.push_back(cnt);
	}
	int code_len = code.size();
	fprintf(fp, "%d", code_len);

	for (int i = 0; i < code_len; i++) {
		fprintf(fp, " %d", code[i]);
	}
	fprintf(fp, "\n");
}

void hirgc::saveOtherData(FILE *fp) { // write auxiliary information to file
	fprintf(fp, "%s\n", meta_data);//identifier
	//------------------------
	writeRunLengthCoding(fp, line_break_len, line_break_vec);//length of shor seuqences
	//------------------------
	fprintf(fp, "%d", pos_vec_len); //intervals of lower-case letters
	for (int i = 0; i < pos_vec_len; i++) {
		fprintf(fp, " %d %d", pos_vec[i].begin, pos_vec[i].length);
	}
	//------------------------
	fprintf(fp, "\n%d", n_vec_len); //intervals of the letter 'N'
	for (int i = 0; i < n_vec_len; i++) {
		fprintf(fp, " %d %d", n_vec[i].begin, n_vec[i].length);
	}
	//------------------------
	fprintf(fp, "\n%d", other_char_len); //other characters
	if (other_char_len > 0) {
		int flag[30];
		for (int i = 0; i < 26; i++) {
			flag[i] = -1;
		}
		vector<int> arr;
	
		for (int i = 0; i < other_char_len; i++) {
			fprintf(fp, " %d", other_char_vec[i].pos);

			int temp = other_char_vec[i].ch;
			if (flag[temp] == -1) {
				arr.push_back(temp);
				flag[temp] = arr.size()-1;
			}
		}

		//save other char information
		int size = arr.size();
		fprintf(fp, " %d", size);
		for (int i = 0; i < size; i++) {
			fprintf(fp, " %d", arr[i]);
		}

		if (size < 10) {
			fprintf(fp, " ");
			for (int i = 0; i < other_char_len; i++) {
				fprintf(fp, "%d", flag[other_char_vec[i].ch]);
			}
		} else {
			int bit_num = ceil(log(size)/log(2));
			int v_num = floor(32.0/bit_num);

			for (int i = 0; i < other_char_len; ) {
				// fprintf(fp, "%d ", flag[other_char_vec[i].ch]);
				unsigned int v = 0;
				for (int j = 0; j < v_num && i < other_char_len; j++, i++) {
					v <<= bit_num;
					v += flag[other_char_vec[i].ch];
				}
				fprintf(fp, " %u", v);
			}
		}
	}
	fprintf(fp, "\n");
}

void hirgc::preProcessRef() { // construction of hash table
	for (int i = 0; i < max_arr_num; i++) {//initial entries
		point[i] = -1;
	}

	// memset(point, -1, sizeof(int)*max_arr_num);

	uint64_t value = 0;
	for (int k = sub_str_num - 1; k >= 0; k--) {
		value <<= 2;
		value += ref_seq_code[k];
	}
	int id = value&(uint64_t)(max_arr_num-1);
	loc[0] = point[id];
	point[id] = 0;

	int step_len = ref_seq_len - sub_str_num + 1;
	int shift_bit_num = (sub_str_num*2-2);
	int one_sub_str = sub_str_num - 1;

	for (int i = 1; i < step_len; i++) {
		value >>= 2;
		value += ((uint64_t)ref_seq_code[i + one_sub_str]<<shift_bit_num);
		
		id = value&(uint64_t)(max_arr_num-1);
		loc[i] = point[id];
		point[id] = i;	
	}	
}

void hirgc::searchMatch(char *refFile) { // greedy matching
	FILE *fp = fopen(refFile, "w");
	if (NULL == fp) {
		printf("ERROR! fail to open file %s\n", refFile);
		return;
	}

	saveOtherData(fp);

	int pre_dismatch_point = 0, pre_pos = 0;
	int step_len = tar_seq_len - sub_str_num + 1;
	int max_length, max_k;

	int dis_str_len = 0, i, id, k, ref_idx, tar_idx, length, cur_pos;;
	uint64_t tar_value;

	for (i = 0; i < step_len; i++) {
		tar_value = 0;
		for (k = sub_str_num - 1; k >= 0; k--) {
			tar_value <<= 2;
			tar_value += tar_seq_code[i+k];
		}

		id = point[tar_value&(uint64_t)(max_arr_num-1)];
		if (id > -1) {
			max_length = -1;
			max_k = -1;
			for (k = id; k != -1; k = loc[k]) {
					//[pos[k], pos[k]+sub_str_num-1] //at least max_arr_num_bit_shift
				ref_idx = k + max_arr_num_bit_shift;//loc[k].pos == k
				tar_idx = i + max_arr_num_bit_shift;
				length = max_arr_num_bit_shift;
				while (ref_idx < ref_seq_len && tar_idx < tar_seq_len && ref_seq_code[ref_idx++] == tar_seq_code[tar_idx++]) {//extend current match
					length++;
				}
				if (length >= sub_str_num && length > max_length) 
				{
					max_length = length;
					max_k = k;
				}
			}
			if (max_k > -1) { //if a match is found
				//first print mismatch substring
				if (dis_str_len > 0) {
					dismatched_str[dis_str_len] = '\0';
					fprintf(fp, "%s\n", dismatched_str);
					dis_str_len = 0;
				}
				//then print match substring
				cur_pos = max_k - pre_pos;
				pre_pos = max_k + max_length;

				fprintf(fp, "%d %d\n", cur_pos, max_length - sub_str_num);

				i += max_length;
				pre_dismatch_point = i;
				if (i < tar_seq_len) {
					dismatched_str[dis_str_len++] = '0'+tar_seq_code[i];
				}
				continue;
			}
		}
		dismatched_str[dis_str_len++] = '0'+tar_seq_code[i];
	}

	for(; i < tar_seq_len; i++) {
		dismatched_str[dis_str_len++] = '0'+tar_seq_code[i];
	}
	if (dis_str_len > 0) {
		dismatched_str[dis_str_len] = '\0';
		fprintf(fp, "%s\n", dismatched_str);
	}
	fclose(fp);
	// gettimeofday(&end,NULL);
	// timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
	// printf("compress timer = %lf ms\n", timer/1000.0);
}

void hirgc::write_time_record(char *res_folder)
{
	char res[1024];
	sprintf(res,"%s/time.txt",res_folder);

	gettimeofday(&t_end,NULL);
	 timer = 1000000 * (t_end.tv_sec - t_start.tv_sec) + t_end.tv_usec - t_start.tv_usec;

	FILE *fp = fopen(res,"a");	
	fprintf(fp,"compress timer = %lf ms; %lf min, thread %d\n", timer/1000.0, timer/1000.0/1000.0/60.0,tid);
	fclose(fp);
}

void hirgc::hirgc_single_ref_compress(char *result_folder, char *ref_file)
{
	gettimeofday(&t_start,NULL);
	
	initial();
	printf("single ref compress\n");
	char ref[1024];
	char tar[1024];
	char res[1024];
	char cmd[1024];

	if(stat( result_folder, &info ) != 0)
	{
		sprintf(cmd,"mkdir %s",result_folder);
		system(cmd);
	}
	
	sprintf(ref,"%s",ref_file);	
	readRefFile(ref);
	preProcessRef();	

	int count = 0;
	for(int i = 0; i < genome_list.size(); i++)
	{
		if ( i >= start && i < end) 
		{
			//printf("single compress, file %d from thread %d\n",i,tid);
			sprintf(tar,"%s",genome_list[i].c_str());
			readTarFile(tar);
			sprintf(res, "%s/%s_ref_%s",result_folder,strip_string(genome_list[i].c_str()).c_str(),strip_string(string(ref_file).c_str()).c_str() );
			// compressFile(ref, tar, res);
			if(stat( res, &info ) == 0)
			{	
				sprintf(cmd,"rm -r %s",res);
				system(cmd);
			}
			searchMatch(res);
		}			
	}

	clear();
	return;	

}

void hirgc::hirgc_cluster_compress(char *result_folder)
{
	gettimeofday(&t_start,NULL);
	int i,j;
	initial();
	int clusters_size[centroid_list.size()];
	char ref[1024];
	char tar[1024];
	char res[1024];
	char cmd[1024];
	
	if(stat( result_folder, &info ) != 0)
	{
		sprintf(cmd,"mkdir %s",result_folder);
		system(cmd);
	}


	int final_ref_id = 0, count, max_count = 0, task_count=0;
	int current_ref_id = -1;

    for(i = 0; i < centroid_list.size(); i++)
	{
		if(task_count >= end)
			break;

		count = 0;
		for(j = 0; j < cluster_list.size(); j++)
		{

			if(cluster_list[j]==i && j != centroid_list[i])
			{

				if( task_count >= start && task_count < end)
				{
					if (current_ref_id != i)
					{
						sprintf(ref,"%s",genome_list[centroid_list[i]].c_str());
						printf("thread %d read ref %s\n",tid,ref);
						readRefFile(ref);
						preProcessRef();						
						current_ref_id = i;
					}


//					printf("thread %d compress file %d %s, cluster %d\n",tid,j, genome_list[j].c_str(),i);
					sprintf(tar,"%s",genome_list[j].c_str());
					readTarFile(tar);
					sprintf(res, "%s/%s_ref_%s", result_folder, strip_string(genome_list[j].c_str()).c_str(),
						strip_string(genome_list[centroid_list[i]].c_str()).c_str());
				// compressFile(ref, tar, res);
				//	if(stat( res, &info ) == 0)	
				//	{
				//		sprintf(cmd,"rm -r %s",res);
				//		system(cmd);
				//	}
					searchMatch(res);	
				}
				task_count++;	
				count++;		
			}
		}

		// check which cluster with the most sampples
		if(count > max_count)
		{
			max_count = count;
			final_ref_id = centroid_list[i];  
		}	
	} 

	if(task_count >= end)
	{
		clear();
		return;
	}
	sprintf(ref,"%s",genome_list[final_ref_id].c_str());
	readRefFile(ref);
	preProcessRef();
    int tar_id;	
	
	for(i = 0; i < centroid_list.size(); i++)
	{

		tar_id = centroid_list[i];
		if( tar_id != final_ref_id)
		{
			if ( task_count >= start && task_count < end )
			{
				printf("compress ref file %d with final ref %d\n", tar_id,final_ref_id);
				sprintf(tar,"%s",genome_list[tar_id].c_str());
				readTarFile(tar);
				sprintf(res, "%s/%s_ref_%s", result_folder, strip_string(genome_list[tar_id].c_str()).c_str(),
						strip_string(genome_list[final_ref_id].c_str()).c_str());
				// compressFile(ref, tar, res);
				//if(stat( res, &info ) == 0)	
				//{
				//	sprintf(cmd,"rm -r %s",res);
				//	system(cmd);
				//}
				searchMatch(res);	
			}
			task_count++;
				
		}

	}
	clear();
	return;	

}



