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
# include <dirent.h>
# include <iterator>
# include <sys/stat.h>
# include <map>
# include <numeric>
# include "Eigen/Eigenvalues" 
# include "Eigen/Dense"
# include "Eigen/Geometry"
# include <assert.h>

# include "dlib-19.8/dlib/clustering.h"
# include "dlib-19.8/dlib/statistics.h"
#define MAX(a,b) ((a) > (b) ? (a) : (b))
using namespace std;
//using namespace dlib;
//using Eigen::MatrixXd;
//using Eigen::VectorXd;
#define MIN(a,b) ((a) < (b) ? (a) : (b))
template<typename A, typename B>
std::pair<B,A> flip_pair(const std::pair<A,B> &p)
{
    return std::pair<B,A>(p.second, p.first);
}
// sort map with accessending order based on value(not key)
template<typename A, typename B>
std::multimap<B,A> flip_map(const std::map<A,B> &src)
{
    std::multimap<B,A> dst;
    std::transform(src.begin(), src.end(), std::inserter(dst, dst.begin()), 
                   flip_pair<A,B>);
    return dst;
}

struct stat info;
class matrix
{

public:
  int col,row;
  double *data;    
  matrix(int m_row,int m_col): row(m_row) , col(m_col) 
  {
    data = new double[col*row];
  }

  double *ith_row(size_t i)
  { return (data+i*col); }

  double operator()(size_t i, size_t j)
  { return data[i*col+j]; }

  void clear()
  { delete data;}

  ~matrix()
  {
    delete data;
  }

};



int *genome_sequence;

std::vector<string> dir_list;
std::vector<string> file_list;

static int agctIndex(char ch) 
{ //encoding rule
    ch = toupper(ch);
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

/**std::vector<string> get_directories(const string& s)
{
    std::vector<string> r;
    for(auto& p : filesystem::recursive_directory_iterator(s))
        if(p.status().type() == filesystem::file_type::directory)
            r.push_back(p.path().string());
    return r;
}**/

int max_index(std::vector<double> v)
{
    double max = v[0];
    int index = 0;
    for( int i = 0; i < v.size(); i++)
    {
        if( v[i] > max)
        {
            max = v[i];
            index = i;
        }
    }

    return index;
}

int max_index(std::vector<double> v, std::vector<int> existing_index)
{
    double max;
    int index = -1;
    for( int i = 0; i < v.size(); i++)
    {
        if ( find(existing_index.begin(),existing_index.end(),i) != existing_index.end() )
            continue;
        if (index == -1)
        {
            index = i;
            max = v[i];
        }
        
        else if( v[i] > max)
        {
            max = v[i];
            index = i;
        }
    }

    return index;
}

double poly_kernel(Eigen::VectorXd x, Eigen::VectorXd y)
{
    int d = 2;
    float c = 1.0;
    double result = 0.0;
    /**for  (int i = 0; i < len; i++)
        result += x[i] * y[i];
    **/
    result  = x.dot(y);
    result += c;
    return pow(result,2);
}

double mountain_function(Eigen::VectorXd x, Eigen::VectorXd y, float alpha = 5.4)
{
    double exponent = 0.0;

    exponent = -poly_kernel(x,x) + 2 * poly_kernel(x,y) -poly_kernel(y,y);
    exponent = exponent * alpha;
    return  exp(exponent); 
}

std::vector<int> subtractive_clustering( Eigen::MatrixXd *input_matrix)
{
    float alpha = 1.0;
    float beta = 1.5;
    std::vector<int> centroids;
    std::vector<double> pontential;
    int row = input_matrix->rows();
    int col = input_matrix->cols();

    Eigen::VectorXd row1,row2;

    for(int i = 0; i < row; i++)
    {
        double mountain_coefficient = 0.0;
        row1 = input_matrix->row(i);
        for(int j = 0; j < row; j++)
        {
            row2 = input_matrix->row(j);
            mountain_coefficient += mountain_function(row1,row2,alpha);
        }
        pontential.push_back(mountain_coefficient);
    }

    int current_centroid = max_index(pontential);
    int max_clusters =   0.01*row > 10 ? 0.01*row : 10;
    double thresshold = 0.3 * pontential[current_centroid];
    centroids.push_back(current_centroid);

    int k = 0.09 * row;
    for(int i = 1; i < k; i++)
    {
        row2 = input_matrix->row(current_centroid);
        for ( int row_index = 0; row_index < row; row_index++ )
        {
            if(find(centroids.begin(),centroids.end(),row_index) == centroids.end())
            {
                row1 = input_matrix->row(row_index);
                pontential[row_index] -= mountain_function(row1,row2,alpha);
            }    
        }
        current_centroid = max_index(pontential,centroids);
        if ( pontential[current_centroid] < thresshold || i > max_clusters )
            break;
    
        centroids.push_back(current_centroid);
    }

    return centroids;
}

void count_tuples(string f_name, Eigen::MatrixXd *m,int row_index, int len, int t_size)
{
    int col = (*m).cols();
    FILE *fp = fopen(f_name.c_str(), "r");

    for(int i = 0; i < col; i++)
    {
        (*m)(row_index,i) = 0.0;
    }
    printf("finish init\n");
    
    char meta_data[1024], chr, ch[1024];

    int letters_len = 0, n_letters_len = 0, index, ch_len, _tar_seq_len = 0;
    int tuple_count,count = 0;

    fgets(meta_data, 1024, fp);
    while (fscanf(fp, "%s", ch) != EOF && count <= len ) {
        ch_len = strlen(ch);

        for (int j = 0; j < ch_len; j++) 
        {
            index = agctIndex(ch[j]);

            if (index != 4) 
                genome_sequence[count++] = index;
        }
                        
    }   
    fclose(fp);
    int tuple_value = 0, value = 0;
    for (int k = t_size - 1; k >= 0; k--) {
        tuple_value <<= 2;
        tuple_value += genome_sequence[k];
    }

    (*m)(row_index,tuple_value) += 1;

    int step_len = len - t_size + 1;
    int shift_bit_num = (t_size*2-2);
    int one_sub_str = t_size - 1;
    for (int i = 1; i < step_len; i++) {
        tuple_value >>= 2;
        tuple_value += genome_sequence[i + one_sub_str]<<shift_bit_num;
        (*m)(row_index,tuple_value) += 1;   
    }    
    return;
}

long double sum_of_variance(Eigen::MatrixXd *input_matrix)
{
    int row = input_matrix->rows();
    int col = input_matrix->cols();
    Eigen::VectorXd new_vector;
    double sum_var = 0.0;

    FILE *fp = fopen("var_record.txt","w");
    for(int i = 0; i < col; i++)
    {
        //current_column.clear();
        new_vector = input_matrix->col(i);
        std::vector<double> current_column(new_vector.data(),new_vector.data() + new_vector.size() ); 

        double dimension = (double)current_column.size();
        double sum = accumulate(current_column.begin(), current_column.end(), 0.0);
        double mean = sum / dimension;

        double sq_sum = std::inner_product(current_column.begin(), current_column.end(), current_column.begin(), 0.0);    
        double variance = sq_sum/dimension - (mean * mean);
        if(variance > 0)
            sum_var = sum_var+variance;
        fprintf(fp,"i: %d, sample size: %lf, sq_sum %lf, mean %lf, variance %lf, sum_var %lf\n"
            ,i, dimension, sq_sum,mean, variance,sum_var);        
    }
    printf("sum of var is %lf\n",sum_var);
    return sum_var;

}

/**
Eigen::MatrixXd PCA(Eigen::MatrixXd *input_matrix, int n_component = 100)
{
    int row = input_matrix->rows();
    int col = input_matrix->cols();

    //Eigen::MatrixXd C = input_matrix->adjoint() * (*input_matrix);
    Eigen::MatrixXd meanval = input_matrix->colwise().mean();
    Eigen::RowVectorXd meanvecRow = meanval;
    input_matrix->rowwise() -= meanvecRow;   
    Eigen::MatrixXd C = input_matrix->transpose() * (*input_matrix);
    printf("dimension of scatter matrix %d %d\n",C.rows(),C.cols());
    for(int i = 0; i < C.rows(); i++ )
    {
        for(int j = 0; j < C.rows(); j++ )
        {
            if( std::isnan( C(i,j) ) )
                C(i,j) = 1.0;
        }
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig;
    //Eigen::EigenSolver<Eigen::MatrixXd> eig;
    printf("compute eigen vectors:\n");
    eig.compute(C);    
    Eigen::MatrixXd vec = eig.eigenvectors().real();
    Eigen::MatrixXd val = eig.eigenvalues().real();

    for( int i = 0; i < vec.rows(); i++)
    {
        for(int j = 0; j < vec.cols(); j++)
        {
            if( std::isnan(vec(i,j)) )
                vec(i,j) = 0.0;
            vec(i,j) = vec(i,j) * 100;
        }
    }    

    int dim = 0;
    for (int i = val.rows()-1; i >= 0; --i)
    {
        dim = i;
        if( val(i,0) <= 0.0001 )
            break;
    }
    dim =  (val.rows()-dim) >= n_component ? n_component : dim;

    Eigen::MatrixXd res = (*input_matrix) * vec.rightCols(val.rows() - dim);
    printf("dimension of reduced matrix: %d * %d\n",res.rows(), res.cols());

    return res;
} **/


Eigen::MatrixXd PCA(Eigen::MatrixXd *input_matrix, int n_component = 100)
{
    Eigen::MatrixXd result;
    int row = input_matrix->rows();
    int col = input_matrix->cols();

    FILE *fp = fopen("selected_features.txt","w");
    for(int i = 0; i < row ; i++)
    {
        for(int j = 0; j < col; j++)
            fprintf(fp, "%lf ", (*input_matrix)(i,j) );
        fprintf(fp,"\n");
    }
    fclose(fp);

    Eigen::VectorXd new_vector;
    Eigen::VectorXd col_mean(col);

    Eigen::MatrixXd center_matrix(*input_matrix);
    for(int i = 0; i <col; i++)
    {
        new_vector = input_matrix->col(i);
        std::vector<double> current_column(new_vector.data(),new_vector.data() + new_vector.size() ); 
        double sum = accumulate(current_column.begin(), current_column.end(), 0.0);
        double mean = sum / current_column.size();
        col_mean(i) = mean; 
        for(int j = 0; j < row; j++)
            center_matrix(j,i) -= mean;
    }

    Eigen::MatrixXd scatter_matrix = center_matrix.transpose() * center_matrix;
    for(int i = 0; i < col; i++)
    {
        for(int j = 0; j < col ; j++)
        {
            if( std::isnan( scatter_matrix(i,j) ))
                scatter_matrix(i,j) = 0.0;
        }
    }
    printf("dimension of scatter matrix:%d %d\n",scatter_matrix.rows(),scatter_matrix.cols());
    Eigen::EigenSolver<Eigen::MatrixXd> es(scatter_matrix); 

    //the princpal component of
    Eigen::VectorXd evalue_real = es.eigenvalues().real();
    Eigen::MatrixXd evector_real = es.eigenvectors().real();
    printf("egien matrix solved\n");
    int count = -1;
    for(int i = 0; i < col; i++)
    {
        if( evalue_real(i) < 0.0001 )
        {
            count = i;
            break;
        }
    }

    if( count != -1 && count < n_component)
    {
        printf("since the %dth egien value is negative\n",count);
        n_component = count;
    }
    Eigen::MatrixXd PC;

    printf("square root of eigen values:\n");

    PC = es.eigenvectors().real().leftCols(n_component);
    result = (*input_matrix) * PC;
    //result = PC.transpose() * (input_matrix->transpose() );
    //result = result.transpose();
    printf("dimension of result matrix:%d %d\n",result.rows(),result.cols());

    return result;
} 

/**
Eigen::MatrixXd PCA(Eigen::MatrixXd *input_matrix, int n_component = 100)
{
    Eigen::MatrixXd result;
    int row = input_matrix->rows();
    int col = input_matrix->cols();

    FILE *fp = fopen("selected_features.txt","w");
    for(int i = 0; i < row ; i++)
    {
        for(int j = 0; j < col; j++)
            fprintf(fp, "%lf ", (*input_matrix)(i,j) );
        fprintf(fp,"\n");
    }
    fclose(fp);

    Eigen::VectorXd new_vector;
    Eigen::VectorXd col_mean(col);

    Eigen::MatrixXd center_matrix(*input_matrix);
    for(int i = 0; i <col; i++)
    {
        new_vector = input_matrix->col(i);
        std::vector<double> current_column(new_vector.data(),new_vector.data() + new_vector.size() ); 
        double sum = accumulate(current_column.begin(), current_column.end(), 0.0);
        double mean = sum / current_column.size();
        col_mean(i) = mean; 
        for(int j = 0; j < row; j++)
            center_matrix(j,i) -= mean;
    }

    Eigen::MatrixXd scatter_matrix = center_matrix.transpose() * center_matrix;
    printf("dimension of scatter matrix:%d %d\n",scatter_matrix.rows(),scatter_matrix.cols());
    Eigen::EigenSolver<Eigen::MatrixXd> es(scatter_matrix); 


    //the princpal component of A

    Eigen::VectorXd evalue_real = es.eigenvalues().real();
    Eigen::MatrixXd evector_real = es.eigenvectors().real();
    printf("egien matrix solved\n");
    int count = -1;
    for(int i = 0; i < col; i++)
    {
        if( evalue_real(i) < 0 )
        {
            count = i;
            break;
        }
    }

    if( count != -1 && count < n_component)
    {
        printf("since the %dth egien value is negative\n",count);
        n_component = count;
    }
    Eigen::MatrixXd PC = Eigen::MatrixXd::Zero(col,n_component);

    printf("square root of eigen values:\n");

    PC = es.eigenvectors().real().leftCols(n_component);
    result = (*input_matrix) * PC;
    //result = PC.transpose() * (input_matrix->transpose() );
    //result = result.transpose();
    printf("dimension of result matrix:%d %d\n",result.rows(),result.cols());

    return result;
} **/

Eigen::MatrixXd similarity(Eigen::MatrixXd *input_matrix)
{
    double sum_var = sum_of_variance(input_matrix);
    double factor = 0.1;
    sum_var *= factor;

    int row = input_matrix->rows();
    int col = input_matrix->cols();
    Eigen::MatrixXd s_matrix(row,row);
    Eigen::VectorXd row1,row2;

    long double s_score;

    FILE *fp = fopen("d_matrix.txt","w");
    FILE *fp2 = fopen("s_matrix.txt","w");

    for(int i = 0; i < row; i++)
    {
        row1 = input_matrix->row(i);
        {
            for(int j = i+1; j < row; j++)
            {
                row2 = input_matrix->row(j);
                double sum = 0.0;
                double buffer;
                for(int k = 0; k < col; k++ )
                {
                    buffer = pow(row1(k) - row2(k),2);
                    if(buffer > 0.0)
                        sum += buffer;
                }
                fprintf(fp,"%lf ",sum);
                //s_score = (long double)( sum/(2*sum_var) ); 
                s_score = sum/(2*sum_var); 
                s_score = exp(-s_score);
                fprintf(fp2,"%Lf ",s_score);
                s_matrix(i,j) = s_score;
                s_matrix(j,i) = s_score;
            }
        }
        fprintf(fp,"\n");
        fprintf(fp2,"\n");
    }
    fclose(fp);
    fclose(fp2);
    return s_matrix;

}



Eigen::MatrixXd SF_selection(Eigen::MatrixXd *input_matrix)
{
    int row = input_matrix->rows();
    int col = input_matrix->cols();    
    int selected_number = MIN(1000,0.2 * col);
    Eigen::MatrixXd selected_features(input_matrix->rows(),selected_number);

    Eigen::MatrixXd w_matrix = similarity(input_matrix);
    if(w_matrix.rows() != w_matrix.cols() )
    {
        printf("warning, weight matrix is not square matrix\n");
        exit(0);
    }

    Eigen::MatrixXd d_matrix = Eigen::MatrixXd::Zero(w_matrix.rows(),w_matrix.cols());
    //power -1/2 of d_matrix
    Eigen::MatrixXd d2_matrix = Eigen::MatrixXd::Zero(w_matrix.rows(),w_matrix.cols());
    for(int i = 0; i < w_matrix.rows(); i++)
    {
        for(int j = 0; j < w_matrix.cols(); j++)
            d_matrix(i,i) += w_matrix(i,j);

        d2_matrix(i,i) = 1/sqrt( d_matrix(i,i) );
  //      cout << "D:" << d_matrix(i,i) << "  D-0.5:" << d2_matrix(i,i) << endl;
    }

    Eigen::MatrixXd l_matrix = d_matrix - w_matrix; 
    printf("dimesnion of L: %d %d\n",l_matrix.rows(),l_matrix.cols());
    Eigen::MatrixXd Lap_matrix = d2_matrix * l_matrix * d2_matrix;

    std::map<int, double> ph1_map;
    Eigen::VectorXd fi; 
    for(int i = 0; i < input_matrix->cols(); i++)
    {
        fi = input_matrix->col(i);
        ph1_map.insert( pair<int,double>(i, ((fi.transpose() * l_matrix * fi ) /  ( fi.transpose() * d_matrix * fi))(0,0) ) );

    }

    std::multimap<double, int> ph1_sort = flip_map(ph1_map);
    std::multimap<double, int>::iterator it = ph1_sort.begin();
    int feature_index;
    for(int i = 0; i < selected_number; i++)
    {
        printf("feature value %lf, index %d\n", it->first,it->second);
        feature_index = it->second;
        selected_features.col(i) = input_matrix->col(feature_index);
        //selected_features.col(i) = input_matrix->col(feature_index);
        ++it;
    }
    exit(0);
    return selected_features;
}


Eigen::MatrixXd frequency_matrix(string src, int len = 1000000,int t_size = 7)
{
    int file_number = file_list.size();
    int frequnecy_dimension = 1 << (2*t_size);
    Eigen::MatrixXd result(file_number,frequnecy_dimension);

    int row = result.rows();
    int col = result.cols();
    char buffer[1024];
    genome_sequence = new int[len+1024];
    for(int row_index = 0; row_index < file_number; row_index++)
    {

        sprintf(buffer,"%s%s",src.c_str(),file_list[row_index].c_str());
        FILE *fp = fopen(buffer, "r");
        for(int i = 0; i < col; i++)
            result(row_index,i) = 0.0;
    
        char meta_data[1024], chr, ch[1024];

        int letters_len = 0, n_letters_len = 0, index, ch_len, _tar_seq_len = 0;
        int tuple_count,count = 0;

        fgets(meta_data, 1024, fp);
    
        while (fscanf(fp, "%s", ch) != EOF && count <= len ) 
        {
            ch_len = strlen(ch);

            for (int j = 0; j < ch_len; j++) 
            {
                index = agctIndex(ch[j]);

                if (index != 4) 
                    genome_sequence[count++] = index;
            }
                        
        }   
        fclose(fp);
        int tuple_value = 0, value = 0;
        for (int k = t_size - 1; k >= 0; k--) {
            tuple_value <<= 2;
            tuple_value += genome_sequence[k];
        }

        result(row_index,tuple_value) += 1;

        int step_len = len - t_size + 1;
        int shift_bit_num = (t_size*2-2);
        int one_sub_str = t_size - 1;
        for (int i = 1; i < step_len; i++) {
            tuple_value >>= 2;
            tuple_value += genome_sequence[i + one_sub_str]<<shift_bit_num;
            result(row_index,tuple_value) += 1;   
        }       
    }
    delete genome_sequence;
    return result;
}

Eigen::MatrixXd normalize(Eigen::MatrixXd *input_matrix)
{
    int row = input_matrix->rows();
    int col = input_matrix->cols();

    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(row,col);
    double max,min, deno;
    for(int i = 0; i < col; i++)
    {
        Eigen::VectorXd current_column = input_matrix->col(i);
        max = current_column.maxCoeff();
        min = current_column.minCoeff();
        deno = max - min;
        if (max == min)
            deno = max;

        for(int j = 0; j < row; j++)
            result(j,i) = (current_column(j) - min)/deno;

    }

    return result;
}


double distance_between_vector(Eigen::VectorXd row1, Eigen::VectorXd row2)
{
    double sum = 0.0;
    for(int k = 0; k < row1.size(); k++ )
        sum += pow(row1(k) - row2(k),2);

    return sum;
}

int minmum_distance_index( std::vector<Eigen::VectorXd> vectors, Eigen::VectorXd row1)
{
    
    double distance = 0.0, min_distance = -1.0;
    int min_index = -1;
    
    Eigen::VectorXd row2;
    for(int i = 0; i < vectors.size(); i++ )
    {
        row2 = vectors[i];
        distance = distance_between_vector(row1,row2);
        if( distance < min_distance || min_distance < 0 )
        {
            min_index = i;
            min_distance = distance;
        }
    }
    return min_index;
}

int minmum_distance_index( Eigen::MatrixXd *input_matrix, std::vector<int> indices,int row_id)
{
    int row = input_matrix->rows();
    int col = input_matrix->cols();

    double distance = 0.0, min_distance = -1.0;
    int min_index = -1;
    
    Eigen::VectorXd row1,row2;
    row1 = input_matrix->row(row_id);
    for(int i = 0; i < indices.size(); i++ )
    {
        row2 = input_matrix->row(i);
        distance = distance_between_vector(row1,row2);
        if( distance < min_distance || min_distance < 0 )
        {
            min_index = i;
            min_distance = distance;
        }
    }
    return min_index;
}

std::vector<Eigen::VectorXd> mean_centroids( Eigen::MatrixXd *input_matrix, std::vector<int> cluster_indices)
{
    int col = input_matrix->cols();
    std::vector<Eigen::VectorXd> vectors;
    std::vector<int> cluster_size;
    int cluster_number = *( max_element(cluster_indices.begin(), cluster_indices.end()) ) + 1;
    int cluster_index;
    
    for(int i = 0; i < cluster_number; i++)
    {
        vectors.push_back(Eigen::VectorXd::Zero(col));
        cluster_size.push_back(0);
    }

    for (int i = 0; i < cluster_indices.size(); i++)
    {
        cluster_index = cluster_indices[i];
        vectors[ cluster_index ] += input_matrix->row(i);
        cluster_size[ cluster_index ] += 1;
    }

    for(int i = 0; i < cluster_number; i++)
    {
        printf("cluster %d ,size %d\n",i,cluster_size[i]);
        vectors[i] /= cluster_size[i];
    }

    return vectors;
}

std::vector<int> centroids_index( Eigen::MatrixXd *input_matrix, std::vector<Eigen::VectorXd> vectors,std::vector<int> cluster_indices)
{
    std::vector<int> centroid_indices;
    std::vector<double> minmum_distance;

    int cluster_number = vectors.size();
    int cluster_index;
    
    for(int i = 0; i < cluster_number; i++)
    {
        centroid_indices.push_back(-1);
        minmum_distance.push_back(-1.0);
    }

    double distance;
    for (int i = 0; i < cluster_indices.size(); i++)
    {
        cluster_index = cluster_indices[i];

        distance = distance_between_vector( vectors[cluster_index],input_matrix->row(i) ); 
        if ( distance < minmum_distance[cluster_index] || minmum_distance[cluster_index] < 0)
        {
            minmum_distance[cluster_index] = distance;
            centroid_indices[cluster_index] = i;
        }
    }

    return centroid_indices;
}


std::tuple<std::vector<int>,std::vector<int>> clustering( Eigen::MatrixXd *input_matrix, std::vector<int> centroids )
{
    std::vector<int> cluster_result;
    int cluster_number = centroids.size();
    printf("number of cluster: %d\n",cluster_number);
    int row = input_matrix->rows();
    int col = input_matrix->cols();

    typedef dlib::matrix<double,0,1> sample_type;
    typedef dlib::radial_basis_kernel<sample_type> kernel_type;

    int max_dictionary = row*0.3;
    dlib::kcentroid<kernel_type> kc(kernel_type(0.1),0.01,max_dictionary);
    dlib::kkmeans<kernel_type> test(kc);
    std::vector<sample_type> samples;
    std::vector<sample_type> initial_centers;
    sample_type m;
    m.set_size(col,1);
    for(int i = 0; i < row; i++)
    {
        for(int j = 0; j < col; j++)
            m(j) = (*input_matrix)(i,j);
        samples.push_back(m);
    }
    test.set_number_of_centers(cluster_number);
    dlib::pick_initial_centers(cluster_number,initial_centers,samples,test.get_kernel());
 /**   for(int i = 0; i < cluster_number; i++)
    {
        int row_index = centroids[i];
        for(int j = 0; j < col; j++)
            m(j) = (*input_matrix)(row_index,j);
        initial_centers.push_back(m);
    } **/


    test.train(samples,initial_centers);
    
    for(int i = 0; i < row; i++)
        cluster_result.push_back( test(samples[i]) );
    
    std::vector<int> result_centroids;
    std::vector<Eigen::VectorXd> centers;
    
    for( int cluster_index = 0; cluster_index < cluster_number; cluster_index++)
    {
        Eigen::VectorXd v = Eigen::VectorXd::Zero(col);
        int count = 0;
        for(int i = 0; i < row; i++)
        {
            if( cluster_result[i] == cluster_index )
            {
                v += input_matrix->row(i);
                count++;
            }
        }
        v /= count;
        centers.push_back(v);
    }

    result_centroids = centroids_index(input_matrix, centers,cluster_result);
    return std::make_tuple(result_centroids,cluster_result);
}

//return centrid index and cluster index
tuple< std::vector<int>,std::vector<int>> kmean_with_initialzation( Eigen::MatrixXd *input_matrix, std::vector<int> centroids)
{
    std::vector<int> cluster_result;
    int row = input_matrix->rows();
    int col = input_matrix->cols();
    int cluster_number = centroids.size();
    printf("number of clusters:%d\n",cluster_number);

    for(int i = 0; i < row; i++)
        cluster_result.push_back(-1);

    for(int i = 0; i < cluster_number; i++)
        cluster_result[ centroids[i] ] = i;

    for(int i = 0; i < row; i++)
    {
        if( cluster_result[i] == -1)
            cluster_result[i] = minmum_distance_index(input_matrix,centroids,i);
    }    


    int file_number = file_list.size();
    FILE *fp = fopen("cluster_result3.txt","w");
    fprintf(fp,"centroid: ");
    char buffer[1024];
    for(int i = 0; i < centroids.size(); i++)
        fprintf(fp,"%d ",centroids[i]);
    fprintf(fp,"\n");
    for(int i = 0; i < file_number; i++)
    {
        sprintf(buffer,"%s",file_list[i].c_str());
        fprintf(fp,"%s %d\n",buffer,cluster_result[i]);
    }
    fclose(fp);
    
    std::vector<int> previous_cluster;
    std::vector<Eigen::VectorXd> centroids_value;

    int iter = 0;
    while( iter < 20 && previous_cluster != cluster_result )
    {
//        printf("iter %d\n",iter);
        previous_cluster = cluster_result;
        centroids_value = mean_centroids(input_matrix,cluster_result);
        for( int i = 0; i < row; i++)
            cluster_result[i] = minmum_distance_index(centroids_value,input_matrix->row(i));
        iter++;
    }

    std::vector<int> centroid_indices = centroids_index(input_matrix,centroids_value,cluster_result);
    return make_tuple(centroid_indices,cluster_result);
}

//length of checked genome_sequence, size of each tuple, y dimesnion of resultant matrix = 4 ^ tuple_size
Eigen::MatrixXd count_frequency(string src, string result_path,int len = 1000000,int tuple_size = 7, int selected_number = 0)
{
    printf("check frequency start\n");
    string default_name = "fa";

    string file_name;
    int file_number = file_list.size();
    int frequnecy_dimension = 1 << (2*tuple_size);

    printf("begin to count tuples\n");
    Eigen::MatrixXd freq = frequency_matrix(src, len, tuple_size);
    if (selected_number == 0)
        selected_number = frequnecy_dimension;

    char buffer[1024];
    genome_sequence = new int[len+1024];


    FILE *fp = fopen(result_path.c_str(),"w");    
    if (NULL == fp) {
        printf("fail to open file %s\n", result_path.c_str());
        exit(0);
    }
    
    std::vector<int> order;
    for(int i = 0; i < freq.cols(); i++)
        order.push_back(i);
  //  random_shuffle(order.begin(),order.end());
    for(int i = 0; i < freq.rows(); i++)
    {
        fprintf(fp,"%s ",file_list[i].c_str());
        for(int j =0; j < selected_number; j++)
            fprintf(fp,"%d ", (int)freq(i,order[j]) );

        fprintf(fp,"\n");
    }

    fclose(fp);
    delete genome_sequence; 
    return freq;
}

int isfa(string name)
{
    if (name.find(".fa") == string::npos)
        return 0;

    if( name.find(".fa") == ( name.length()-strlen(".fa")) || name.find(".fasta") == ( name.length()-strlen(".fasta")) )    
        return 1;
    
    return 0;
}

void search_files(string src)
{
    string file_name;
    int count = 0;
    DIR *dir = opendir(src.c_str());
    struct dirent *file;
    while ( (file = readdir(dir)) != NULL)
    {
        file_name = file->d_name;
        if( isfa(file_name) )
        {
            file_list.push_back(file_name);
//            printf("found file: %s\n",file_name.c_str());
            count++;   
        }
    }
//    printf("finish search\n");
    return;
}

void search_files(string src, string src_list)
{
 //   printf("notice: src list is supposed to be each line contain name of a file\n");
    FILE *src_file = fopen(src_list.c_str(),"r");
    char file_name[1024];
    char final_fp[1024];
    printf("open file %s\n",src_list.c_str());
    while (fscanf(src_file, "%s", file_name) != EOF) 
    {

        sprintf(final_fp,"%s%s",src.c_str(),file_name);
        if(stat( final_fp, &info ) == 0)
        {
            
            //printf("find file %s\n",final_fp);
            file_list.push_back(file_name);
        }
        else { printf("fail to find file:%s\n",final_fp); }
    }
    return;
}

class preprocessor
{
public:
    string src;
    string result_path;
    string src_list;
    int tuple_size, len,kept_dimension;
    std::vector<int> centroid_result;    
    std::vector<int> cluster_result;
    preprocessor(string m_src, string m_result_path = "",string m_src_list = "",int m_len = 1000000,int m_tuple_size = 7,int m_kept_dimension = 0 )
                : src(m_src), result_path(m_result_path), src_list(m_src_list), tuple_size(m_tuple_size), len(m_len), kept_dimension(m_kept_dimension) {}

    void preprocess()
    {

        search_files(src,src_list); 
        printf("read length is set as %d (default 1000000), tuple size if set as %d(default 7)\n",len,tuple_size);
        Eigen::MatrixXd freq = count_frequency(src,result_path,len, tuple_size,kept_dimension);
        Eigen::MatrixXd normalized_data = normalize(&freq);
        Eigen::MatrixXd selected_features = SF_selection(&normalized_data);
        for( int i = 0; i < selected_features.rows(); i++)
        {
            for(int j = 0; j < selected_features.cols(); j++)
            {
                if( std::isnan(selected_features(i,j)) )
                    selected_features(i,j) = 0.0;
            }
        }
        Eigen::MatrixXd data = PCA(&selected_features);
        std::vector<int> centroids = subtractive_clustering(&data);
    
        tie(centroid_result,cluster_result) = clustering(&data,centroids);
        return;
    }

    std::vector<int> get_centroid()
    {
        return centroid_result;
    }

    std::vector<int> get_cluster()
    {
        return cluster_result;
    }
    
};
