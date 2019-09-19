#! /bin/bash
# pom2k for auto-parallel job
#PBS -q s
#PBS -l nodes=1:ppn=20
#PBS -N POM2K
#PBS -o std_out
#PBS -e std_err

 cd $PBS_O_WORKDIR
 export OMP_NUM_THREADS=20

 rm -rf pom2k.exe
 rm -rf error_check.txt
 rm -rf std_err
 rm -rf std_out
#rm -rf case 

 result_dir="/work1/gaoj/POM/pom_seamount/result/" 
 file_dir="$result_dir`date +%Y%m%d%H%M`" 
#if [ ! -x "$result_dir" ];then 
#mkdir $result_dir 
#fi 

 echo $file_dir
 #cd $result_dir
 cd $result_dir
 cd ../main_prog
 echo $file_dir> case.txt

 mkdir $file_dir

 cp -rf /work1/gaoj/POM/pom_seamount/main_prog $file_dir

 cd $file_dir
 cd main_prog
# pwd
 rm -rf pom2k.exe
 ifort -O3 -mcmodel=large -parallel -o pom2k.exe pom2k.f sub*.f 
#-I/ap/netcdf/3.6.3/include -L/ap/netcdf/3.6.3/lib -lnetcdf 

#ifort -parallel -par-report2  pom2k.f sub*.f -I/ap/netcdf/3.6.3/include/ -L/ap/netcdf/3.6.3/lib -lnetcdf -O2

  ./pom2k.exe
#./a.out
  rm -rf fort.*
 cd ../
 mkdir output_pic
