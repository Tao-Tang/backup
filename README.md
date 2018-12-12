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
Process Mode: compute reference for each file in dataset, result is a text file with multiple lines, each line contains two file name, the first is reference sequence and second is target sequence.
```
./ECC p -s <src> -r <result>
```
Compress Mode: Execute reference selection process and compress files based on selection result with HIRGC, output includes <result>.tar <result>_src.txt <result>_ref  <br />
```
./ECC c -s <src> -r <result> -t <thread_number> 
```
Decompress Mode: decompress file
```
./ECC d  <result>
```
```
<src>:  name of src file, each line records the file path of one sequence  
<result>: name of result file  
<thread_number>: number of threads used
```
