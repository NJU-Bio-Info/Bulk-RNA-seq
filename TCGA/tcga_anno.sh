#使用步骤#
#第一步：cd sh所在目录#
#第二步：用chmod和把sh加载进去#
#第三步：用 sed -i 's/\r//' 你的文件.sh 防止因为编写系统不同造成的bug#
#第四步：【仍在这个目录下】输入sh tcga_anno.sh -d <dir> -a <anno>#
#!/bin/bash
help(){
	cat <<EOF
Description
	
	This script is used to annotation gene expression profiling data downloaded from TCGA database.

	A required file: the sample sheet downloaded from TCGA database.
		
	All directory containing RPKM gz file must be placed in a clean directory.

Usage
	tcga_anno.sh -d <dir> -a <anno>#使用的时候要加sh命令在前面#
	
Options
	-d		the directory containing data
	-a		the annotation sample information sheet file

Last updated: 2022/3/7
Owner: ShuiKM

EOF
}

if [ $# = 0 ] || [ $1 = '-h' ] || [ $1 = '--help' ]
then
	help
	exit 1
fi
while getopts "d:a:" OPTION
do
	case $OPTION in
		d) dir=$OPTARG;;#相当于把输进去的d，a分别装到dir和anno变量里#
		a) anno=$OPTARG;;
	esac
done
cd $dir
ls $dir | while read id#这个循环相当于把所有的gz压缩包从文件夹里移出来了#
do
	cd $id
	mv *.gz ../
  cd ../
done
ls $dir | grep -v ".gz$" | while read id
do
	rm -rf $id#这里相当于把移出gz后空的文件夹删掉了#
done
i=1
  ls $dir | while read id#相当于做一个循环把里面的gz文件都读完#
do
	temp=$(cat $anno | awk -v name=$id 'BEGIN{FS="/t"}$2==name{print$8}')
  sample=$(echo $temp | sed 's/ /_/g')
	mv $id $sample"_"$i.txt.gz
	i=$[ $i+1 ]
done

gunzip *.gz

unset temp
unset sample



