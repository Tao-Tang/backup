# ECC
### ECC - What is it?
ECC is a reference selction algorithm for genome sequence set ( normally for a set of fasta files ). 
 <br />

### Compile
download this repository <br />
```
cd ECC
chmod +x make
./make
```
### Usage
Process Mode: compute reference for each file in dataset
```
./ECC process -s <src> -r <result> -t <thread_number> 
```


### Example 
```
mpirun -np 32 2dring -N 1024 -M 1000 -output correaltion_matrix.txt -src sample_data.txt -x 2
```
or
```
srun -n 16 -d 24 2dring -N 1024 -M 1000 -output correaltion_matrix.txt -src sample_data.txt -x 2
```
<br />
