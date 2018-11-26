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
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "hirgc3.cpp"

using namespace std;

class Compressor
{
	protected:
		vector<string> file_list;
		vector<int> cluster_index;
		vector<int> centroids;

	public:
		int compress_count = 0;		
		string result_name;
		string ref_file;
		int thread_number = 1;
		struct stat info;	

		Compressor(string m_r, vector<string> m_file_list, vector<int> m_centroids, vector<int> m_cluster_index):
		result_name(m_r), file_list(m_file_list), centroids(m_centroids), cluster_index(m_cluster_index)
		{	}

		Compressor(string m_r, vector<string> m_file_list, string m_ref):
		result_name(m_r), file_list(m_file_list), ref_file(m_ref)
		{	}		
		void set_thread_number(int m_thread_number) { thread_number = m_thread_number;}
		virtual void compress();
		virtual void cluster_compress();
		virtual void ref_compress();
	    
		virtual ~Compressor() {}

};

void Compressor::compress()
{
	char cmd[1024];
	sprintf(cmd,"rm -rf %s",result_name.c_str() );
	system(cmd);
	sprintf(cmd,"mkdir -p %s",result_name.c_str());
	system(cmd);

	if( cluster_index.size() > 0 )
		cluster_compress();

	else if ( ref_file != "" )
	{
		printf("compress with first file as reference\n");
		ref_compress();
	}

}


void Compressor::cluster_compress()
{
	printf("number of files: %d\n",file_list.size());

	char res_folder[1024];
	sprintf(res_folder,"%s",result_name.c_str());

	thread *threads[thread_number];
	hirgc *hirgcs[thread_number];
	int start = 0, task_number = file_list.size();
	int end = start;
	int avg_task = task_number/thread_number;
	int moduler = task_number%thread_number;

	for (int i = 0; i < thread_number; i++)
	{
		end +=  i < moduler ? avg_task+1 : avg_task;
		hirgcs[i] = new hirgc(file_list,start,end,i,cluster_index,centroids);
		threads[i] = new thread(&hirgc::hirgc_cluster_compress,hirgcs[i],res_folder);
		start = end;
	}

	for (int i = 0; i < thread_number; i++)
	{
		threads[i]->join();
		delete hirgcs[i];
		delete threads[i];
	}

	string ref_name;
	char cmd[1024];
	char current_directory[1024];
	char folder_buffer[1024];
	getcwd(current_directory, 1024);
	chdir(res_folder);		

	for (int i = 0; i < centroids.size(); i++)
	{
		ref_name = strip_string(file_list[centroids[i]]);
		sprintf(folder_buffer,"%s_ref",ref_name.c_str());
		sprintf(cmd,"mkdir %s",folder_buffer);
		system(cmd);                
		sprintf(cmd,"mv *_ref_%s %s",ref_name.c_str(),folder_buffer);
		system(cmd);
		sprintf(cmd,"tar -cf %s.tar -C %s .",folder_buffer,folder_buffer);
		system(cmd);
		sprintf(cmd,"%s/bsc e %s.tar %s.tar.bsc -b64p",current_directory,folder_buffer,folder_buffer);
		printf("bsc command is %s\n",cmd);
		system(cmd);		
	}

	sprintf(cmd, "mkdir bsc");
	system(cmd);
	sprintf(cmd,"mv *.bsc bsc");
	system(cmd);
	sprintf(cmd,"tar -cf %s.tar -C bsc .",res_folder);
	system(cmd);
	chdir(current_directory);
	sprintf(cmd,"rm -rf %s.tar",res_folder);
	system(cmd);	
	sprintf(cmd,"mv %s/%s.tar ./",res_folder,res_folder);
	system(cmd);
}

void Compressor::ref_compress()
{
	printf("number of files: %d\n",file_list.size());
	char cmd[1024];;
	char ref[1024];
	sprintf(ref,"%s",ref_file.c_str());

	char res_folder[1024];
	sprintf(res_folder,"%s_ref",ref);

	sprintf(cmd,"mkdir  %s",res_folder);
	system(cmd);
	
	thread *threads[thread_number];
	hirgc *hirgcs[thread_number];
	int start = 0, task_number = file_list.size();
	int end = start;
	int avg_task = task_number/thread_number;
	int moduler = task_number%thread_number;

	for (int i = 0; i < thread_number; i++)
	{
		end +=  i < moduler ? avg_task+1 : avg_task;
		hirgcs[i] = new hirgc(file_list,start,end,i);
		threads[i] = new thread(&hirgc::hirgc_single_ref_compress,hirgcs[i],res_folder,ref);
		start = end;
	}

	for (int i = 0; i < thread_number; i++)
	{
		threads[i]->join();
		delete hirgcs[i];
		delete threads[i];
	}


	char folder_buffer[1024];

	sprintf(cmd,"tar -cf %s.tar -C %s .",res_folder,res_folder);
	system(cmd);

	sprintf(cmd,"./bsc e %s.tar %s.tar.bsc -b64p",res_folder,res_folder);
	system(cmd);

	sprintf(cmd,"mkdir %s",result_name.c_str());
	system(cmd);
	sprintf(cmd,"mv %s.tar.bsc %s",res_folder,result_name.c_str());
	system(cmd);
	sprintf(cmd,"tar -cf %s.tar -C %s .",result_name.c_str(),result_name.c_str());		
	system(cmd);
	return;
} 