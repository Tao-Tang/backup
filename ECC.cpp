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
# include "hirgc-decompressor.cpp"
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

void generate_src_file(char *result,std::vector<string> file_list)
{
	char file_name[1024];
	sprintf(file_name,"%s_src.txt", result);
	FILE *fp = fopen(file_name,"w");
	for(int i = 1; i < file_list.size(); i++)
	{
		fprintf(fp, "%s %s\n", file_list[0].c_str(), file_list[i].c_str());
	}

	fclose(fp);
	return;
}

string generate_src_file(char *result, std::vector<int> centroid_index, std::vector<int> cluster_result, std::vector<string> file_list)
{
	char file_name[1024];
	printf("number of cluster: %d\n", centroid_index.size());
	sprintf(file_name,"%s_src.txt", result);
	FILE *fp = fopen(file_name,"w");	
	int count, max_count = -1, final_ref_id = -1;
    for(int i = 0; i < centroid_index.size(); i++)
	{
		if ( centroid_index[i] == -1)
			continue;

		count = 0;
		for(int j = 0; j < cluster_result.size(); j++)
		{
			if(cluster_result[j]==i && j != centroid_index[i])
				count++;
		}

		if(count > max_count)
		{
			max_count = count;
			final_ref_id = centroid_index[i];  
		}			
	}
	for(int i = 0; i < centroid_index.size(); i++)
	{
		if ( centroid_index[i] == -1)
			continue;
					
		int tar_id = centroid_index[i];
		if( tar_id != final_ref_id)
			fprintf(fp,"%s %s\n", file_list[ final_ref_id ].c_str(), file_list[ tar_id ].c_str() );
	}	

    for(int i = 0; i < centroid_index.size(); i++)
	{
		if ( centroid_index[i] == -1)
			continue;

		count = 0;
		for(int j = 0; j < cluster_result.size(); j++)
		{
			if(cluster_result[j]==i && j != centroid_index[i])
				fprintf(fp,"%s %s\n",file_list[ centroid_index[i] ].c_str(), file_list[j].c_str() );
		}
	}	
	fclose(fp);
	return file_list[ final_ref_id ] ;
}


void decompress_mode(int argc, char *argv[])
{
	printf("decompress file\n");
	string src_name(argv[2]);
	printf("archive name: %s\n",src_name.c_str());
	char cmd[1024];
	sprintf(cmd, "tar -xvf %s.tar",src_name.c_str());
	system(cmd);
	printf("begin read, src is %s\n",src_name.c_str());
	sprintf(cmd,"%s_src.txt",src_name.c_str());
	FILE *fp = fopen( cmd,"r");
	if( fp == NULL)
	{
		printf("failed to open file\n");
		exit(0);
	}

	std::vector<string> ref_list;
	std::vector<string> tar_list;
	std::vector<string> refs;

	char temp_ch[1024];
	string temp_ref=" ";
	printf("begin to read src\n");
	while(fscanf(fp,"%s",temp_ch) != EOF)
	{
		ref_list.push_back( string(temp_ch) );
		if ( temp_ref != string(temp_ch))
		{
			temp_ref = string(temp_ch);
			refs.push_back(temp_ref);
			printf("add ref %s\n",temp_ref.c_str());
		}
		fscanf(fp,"%s",temp_ch);
		tar_list.push_back( string(temp_ch) );
	}

	for (int i = 0; i < refs.size(); i++)
	{
		sprintf(cmd,"./bsc d  %s_ref.tar.bsc %s.tar",refs[i].c_str(),refs[i].c_str());
		system(cmd);
		sprintf(cmd,"tar -xvf %s.tar",refs[i].c_str());
		system(cmd);
		sprintf(cmd,"rm -xvf %s.tar",refs[i].c_str());
		system(cmd);
	}

	sprintf(cmd,"./bsc d %s_ref %s",src_name.c_str(),ref_list[0].c_str());
	system(cmd);
	deCompressor decomp(ref_list[0]);
	char buffer1[1024];
	char buffer2[1024];
	for(int i  = 0; i < ref_list.size(); i++)
	{
		sprintf(buffer1,"%s",ref_list[i].c_str());
		sprintf(buffer2,"%s",tar_list[i].c_str());
		decomp.decompress(buffer1, buffer2 );
	}

	return;
}

void compress_mode(int argc, char *argv[], bool compress=true)
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
	char *src_file;
	int thread_number = 1, n_count = 1, tuple_size = 7, kept_dimension = 1000 ;
	bool single_ref = false;
	bool idc = false;
	while( n_count < argc )
	{
		buffer = argv[n_count];
		if( strcmp(buffer,"-s" )== 0 )
		{
			src_file = n_count < (argc-1) ? argv[++n_count] : NULL;
			src_list.push_back(string(buffer));
		}

		else if( strcmp(buffer,"-sl") == 0)
		{
			src_file = n_count < (argc-1) ? argv[++n_count] : NULL;
			//src_list = read_list(buffer);
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

		else if (strcmp(buffer,"-tz") == 0)
		{
			buffer = n_count < (argc-1) ? argv[++n_count] : NULL;			
			tuple_size = atoi(buffer);		
		}

		else if (strcmp(buffer,"-kd") == 0)
		{
			buffer = n_count < (argc-1) ? argv[++n_count] : NULL;			
			kept_dimension = atoi(buffer);		
		}
		n_count++;
	}	

	char result_name[1024];
	int len = 1000000; 

	printf("src file is %s\n",src_file);
	std::vector<string> file_list = read_list(string(src_file));

	sprintf(result_name,"%s",result);
	char cmd[1024];
	sprintf(cmd,"rm -rf %s %s.tar",result,result);
	system(cmd);

	if( !single_ref )
	{
		preprocessor processor("","process_record.txt",src_file,len,tuple_size,kept_dimension);
		system("rm -rf process_record.txt");
   		processor.preprocess(); 
		std::vector<int> centroid_index = processor.get_centroid();
		std::vector<int> cluster_result = processor.get_cluster();
	
		string f_ref = generate_src_file(result,centroid_index,cluster_result,file_list);
		if( compress == false )
		{ 
			return;
		}
		sprintf(cmd,"./bsc e %s %s_ref -b64p",f_ref.c_str(),result);
		system(cmd);
		Compressor compressor(result_name, file_list, centroid_index, cluster_result);
		compressor.set_thread_number(thread_number);
		compressor.compress();
	}

	else
	{
		string ref = string(file_list[0]);
		generate_src_file(result_name,file_list);
		file_list.erase( file_list.begin() );
		Compressor compressor(result_name, file_list, ref);
		compressor.set_thread_number(thread_number);
		compressor.compress();
	}
	sprintf(cmd,"rm -rf %s",result);
	system(cmd);
	gettimeofday(&end,NULL);
	timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
	printf("compression takes = %lf ms; %lf min, thread_number = %d\n", timer/1000.0, timer/1000.0/1000.0/60.0,thread_number);		

}





static void show_usage()
{
	printf("v1.0 \n"
		    "./ECC p for reference-target pair selectoin\n"
			"./ECC c for reference-target pair selectoin and compress data via hirgc\n"
			"./ECC d for decompression\n");
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

	string bsc = "bsc";
	ifstream f(bsc.c_str());
	if (!f.good())
	{
		printf("bsc not found\n");
		return 0;
	}

	if ( strcmp(mode,"p")==0 )
	{
		compress_mode(argc,argv,false);
	} 

	else if ( strcmp(mode,"c")==0 )
	{ 
		compress_mode(argc,argv); 
	}
	
	else if ( strcmp(mode,"d")==0 )
	{
		decompress_mode(argc,argv);
	}



	return 0;
}