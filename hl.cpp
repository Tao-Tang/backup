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
# include "process-src4.cpp"
# include "hirgc-compressor4.cpp"

using namespace std;

std::vector<string> read_list(string file_name)
{
	FILE *fp = fopen(file_name.c_str(),"r");
	std::vector<string> result;
	if( fp == NULL)
	{
		printf("failed to open file %s\n",file_name.c_str());
		exit(0);
	}

	char temp_ch[1024];
	while(fscanf(fp,"%s",temp_ch) != EOF)
	{ 
		result.push_back(string(temp_ch));
	}

	return result;
}

void compress_mode(int argc, char *argv[])
{
	printf("compress mode begin:\n");
	struct  timeval  start;
	struct  timeval  end;
	unsigned long timer;
	gettimeofday(&start,NULL);

	if(argc < 3)
	{
		printf("sample command:\n");
		return;
	}

	char *result, *buffer;
	std::vector<string> src_list;
	int thread_number = 1, n_count = 1;
	bool single_ref = false;
	while( n_count < argc )
	{
		buffer = argv[n_count];
		if( strcmp(buffer,"-s" )== 0 )
		{
			buffer = n_count < (argc-1) ? argv[++n_count] : NULL;
			src_list.push_back(string(buffer));
		}

		else if( strcmp(buffer,"-sl") == 0)
		{
			buffer = n_count < (argc-1) ? argv[++n_count] : NULL;
			src_list = read_list(buffer);
		}

		else if( strcmp(buffer,"-r") == 0)
			result = n_count < (argc-1) ? argv[++n_count] : NULL;
		
		else if( strcmp(buffer, "-ref") == 0 )
		{	single_ref = true;   }

		else if( strcmp(buffer,"-t") == 0)
		{
			buffer = n_count < (argc-1) ? argv[++n_count] : NULL;			
			thread_number = atoi(buffer);
		}	

		n_count++;
	}	

	char result_name[1024];
	int len = 1000000, tuple_size = 7, reserved_dimension = 1000; 
	for( int i = 0; i < src_list.size(); i++)
	{
		std::vector<string> file_list = read_list(src_list[i]);

		if( src_list.size() > 1)
		{	sprintf(result_name,"%s/%d",result,i);	}
		else
		{	sprintf(result_name,"%s",result);	}

		if( !single_ref )
		{
			preprocessor processor("","process_record.txt",src_list[i],len,tuple_size,reserved_dimension);
			system("rm -r process_record.txt");
    		processor.preprocess(); 
			std::vector<int> centroid_index = processor.get_centroid();
			std::vector<int> cluster_result = processor.get_cluster();
			printf("list %s ref files: \n", src_list[i].c_str() );
			for( int j = 0; j < centroid_index.size(); j++)
				printf("%s\n",file_list[centroid_index[j] ].c_str() );
		
			Compressor compressor(result_name, file_list, centroid_index, cluster_result);
			compressor.set_thread_number(thread_number);
			compressor.compress();
		}

		else
		{
			string ref = string(file_list[0]);
			file_list.erase( file_list.begin() );
			Compressor compressor(result_name, file_list, ref);
			compressor.set_thread_number(thread_number);
			compressor.compress();
		}

	}
		gettimeofday(&end,NULL);
		timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;

	FILE *fp = fopen("compress_record.txt","a");
	fprintf(fp, "compression takes = %lf ms; %lf min, thread_number = %d\n", timer/1000.0, timer/1000.0/1000.0/60.0,thread_number);		
	fclose(fp);

}


/**
void compress_mode(int argc, char *argv[])
{
	Compressor *compressor;
	printf("enter ./main compress -h to get help\n");
	if(argc < 3)
		return;

	int n_count = 2, thread_number = 1;

	char *folder_list,*cluster_file,*ref_file;
	char *mode = argv[2];


	while( n_count < argc )
	{
		buffer = argv[n_count];

		if( strcmp(buffer,"-h") == 0)
		{
			compress_help();
			return;
		}
		
		else if( strcmp(buffer,"-cluster") == 0)
			compressor = new Cluster_Compressor();
		
		else if( strcmp(buffer,"-single")==0 )
			compressor = new Single_File_Compressor();

		

	}

**/


static void show_usage()
{
	printf("v1.0 \n"
		    "./main process for data pre-processing (genome sequence data in fa or fasta format)\n"
			"./main cluster for clustering based on resultant data of pre-processing\n"
			"./main compress for compression based on clusteirng result via hirgc\n");
}



int main(int argc, char *argv[])
{
//	Py_Initialzie();
//	pName = PyString_FromString("clustering");
//	pModule = PyImport_Import(pName);
	char *mode;

	if( argc > 1)
	{ mode = argv[1]; }

	else 
	{
		show_usage();
		return 0;
	}


	if ( strcmp(mode,"process")==0 )
	{} 
/**	 process_mode(argc,argv); 

	else if ( strcmp(mode,"cluster")==0 )
	{ cluster_mode(argc,argv); } **/

	else if ( strcmp(mode,"compress")==0 )
	{ compress_mode(argc,argv); }
	
	else 
	{
		// printf("arguments error...\n");
		show_usage();
		return 0;
	}
	return 0;
}