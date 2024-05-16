#include <iostream>
#include <stdio.h>
#include <fstream>
#include <cstdlib>
#include <cstring>
using namespace std;
#include "whash.h"
// 基础转换数组，用于将字符（碱基）转换为数字。'A'=0, 'C'=1, 'G'=2, 'T'=3, 非碱基字符转换为4。
uint8_t trans[128]= {
		4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
		4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
		4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
		4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
		4,0,4,1,4,4,4,2,4,4,4,4,4,4,2,4,
		4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,
		4,0,4,1,4,4,4,2,4,4,4,4,4,4,2,4,
		4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4
	};
// 比较两个种子（seed）结构，用于快速排序（qsort）函数。
int comparator(const void *p, const void *q)
{
	seed 	*t = (seed *)p;
	seed 	*f = (seed *)q;
	if (t->sed_value > f->sed_value)
		return 1;
	else {
		if(t->sed_value == f->sed_value) {
			if (t->sed_bkt_pos > f->sed_bkt_pos) {
				return 1;
			} else {
				if(t->sed_bkt_pos == f->sed_bkt_pos)
					return 0;
				else
					return -1;
			}

		}
		else
			return -1;
	}
}
// 将一段基因组序列（碱基字符串）转换为一个数字，用于哈希。
int transfer(char *genome,uint32_t *l2r,uint32_t len_sed)
{

	*l2r = 0;
	uint32_t temp;
	for(uint32_t i = 0;i<len_sed;i++) {
		temp = trans[genome[i]];
		*l2r = *l2r << 2;
		*l2r = *l2r | temp;
	}
	return 0;
}
// 统计去重后的种子数量。
uint32_t statsed(seed *repsed,uint32_t len)
{	
	uint32_t count = 1;
	uint32_t temp_sed_value = repsed[0].sed_value;
	uint32_t temp_sed_bkt_pos = repsed[0].sed_bkt_pos;
	for (uint32_t i=1;i<len;++i) {
		if (repsed[i].sed_value == temp_sed_value && repsed[i].sed_bkt_pos == temp_sed_bkt_pos )
			continue;
		else {
			++count;
			temp_sed_value = repsed[i].sed_value;
			temp_sed_bkt_pos = repsed[i].sed_bkt_pos;			
		}
	}
	return count;

}
// 定义Hash类的成员函数write_hashfile，该函数用于写入哈希文件
// 参数path是哈希文件保存的目录路径，参数genome是待索引的基因组数据，参数len_genome是基因组数据的长度，参数len_sed是哈希种子的长度（通常是k-mer的长度）
int Hash::write_hashfile(char *path,char *genome,uint32_t len_genome, uint32_t len_sed)
{   // 初始化文件流数组，用于保存NUM_FILE个哈希文件
	fstream fp[ NUM_FILE ];
	uint32_t count_line[NUM_FILE]; // 记录每个哈希文件中行的数量
	uint32_t N_count; // 'N'字符计数器，用于处理基因组中的未知碱基
	uint32_t count_real_line; // 实际有效行数（用不到，但代码中有声明）
	char num[4] = {0}; // 数字字符串，用于构造文件名
	char filepath[256] = {0}; // 完整的文件路径
	// 遍历NUM_FILE个文件，为每个文件创建文件流
	for(int i=0;i<NUM_FILE;i++) { // 构造每个文件的路径名
		snprintf(num,4,"%d",i);
		strcpy(filepath,path);
		strcat(filepath,"/");
		strcat(filepath,num);
		//cout<<filepath<<endl; 打开或创建文件，用于读写，如果文件已存在则清空
		fp[i].open(filepath,ios::in|ios::out|ios::trunc|ios::binary);
		if(fp[i].fail()) { // 文件打开失败，打印错误信息并退出
			fprintf(stderr,"Failed to create hash file, now exit");
			exit(1);
		} //else 
			//fprintf(stdout,"open %d successfully")
		count_line[i] = 0; // 初始化行计数器为0
	}
	// 定义掩码，用于处理哈希值和处理'N'的情况
	uint32_t num_pre;
	int move_times = (len_sed - 4)<<1;
	uint32_t mask = 0xffffffff >> (32 - move_times);
	uint32_t mask2 = 0xffffffff >> (32 - (len_sed<<1));
	uint32_t file_seq = 0;
	//change reference start position from none 'N' string
	// 使用transfer函数将基因组的第一个len_sed长度的序列转换成数字，并根据该数字决定文件序号
	transfer(genome,&num_pre,len_sed);
	file_seq = num_pre>>move_times;
	//fp[file_seq]<< (num_pre&mask)<<" "<<0<<endl;
	// 写入第一个哈希值和位置
	uint32_t _value = num_pre&mask;
	uint32_t _pos = 0;

	fp[file_seq].write((char *)&_value,sizeof(uint32_t));
	fp[file_seq].write((char *)&_pos,sizeof(uint32_t));
	// 更新文件中行的数量
	++count_line[file_seq];

	N_count = 0; // 初始化'N'计数器为0
	//fprintf(stderr,"%d\n",len_genome);
	// 遍历基因组序列，将每个子序列哈希到相应的文件中
	for( uint32_t i = 1; i <= len_genome - len_sed; i++ ) {
		//first if genome[i+len - 1 ] == 'n'
		//Yes:N_count ++
		//judge if N_count > N_LIMIT
		// No:   i += len_sed 
		// for j = 0 j+ i = 'n'; ++i);
		// 
		// No: N_count = 0
		// 处理含有'N'的序列
		if ('N' == genome[i + len_sed - 1]) {
			++N_count;
			if (N_count > N_LIMIT) {
				i += len_sed - 1 ;
				uint32_t j; // 寻找下一个非'N'的位置
				for (j=1;j+i<=len_genome-len_sed&& 'N' == genome[j+i];++j);
				//fprintf(stderr,"%d\t",i+j);
				if (i+j > len_genome - len_sed) break; // 如果到达末尾，则终止循环
				i += j;	
				N_count = 0;//here may be wrong cause the i may be len_genome
				transfer(genome + i,&num_pre,len_sed);//don't care if there is 'N' in genome + i + 1 now 为下一段序列计算哈希值
				file_seq = num_pre>>move_times;
				//fp[file_seq]<<(num_pre&mask)<<" "<<(i>>10)<<endl;
				_value = num_pre&mask; // 写入哈希值和位置
				_pos = i>>10;
				
				//cout<<i<<"\t"<<file_seq<<endl;

				fp[file_seq].write((char *)&_value,sizeof(uint32_t)); // 将变量 _value 的内容写入到文件 fp[file_seq] 中
				fp[file_seq].write((char *)&_pos,sizeof(uint32_t)); // 将变量 _pos 的内容写入到文件 fp[file_seq] 中
				
				++count_line[file_seq];	
				continue;
			} 
		} else {
			N_count = 0; // 写入哈希值和位置
		} // 为当前序列计算哈希值，并写入对应的文件
		num_pre = (num_pre<<2&mask2)|trans[genome[i + len_sed - 1]];//it will be corrupted by a problem since num_pre doesn't exist any more 
		file_seq = num_pre>>move_times;
		//fp[file_seq]<<(num_pre&mask)<<" "<<(i>>10)<<endl;
		_value = num_pre&mask;
		_pos = i>>10;
		//cout<<i<<"\t"<<file_seq<<endl;
		fp[file_seq].write((char *)&_value,sizeof(uint32_t));
		fp[file_seq].write((char *)&_pos,sizeof(uint32_t));
		
		++count_line[file_seq];	
	}
	//fprintf(stderr,"creat files successfully\n");
	//
	
	uint32_t amount = 1 << (len_sed<<1);
	// 为哈希表分配内存并初始化
	Hashtab * hashtab = new Hashtab;//may be something wrong;
	if( NULL == hashtab) {
		fprintf(stderr,"Failed when applying for new space");
		exit(1);
	}

	hashtab->pointer = new uint32_t[amount + 1];//may be something wrong;

	if ( NULL == hashtab->pointer ) {
		//cout<<"<<endl;
		fprintf(stderr,"Failed when applying for new space");
		exit(1);

	}

	for(uint32_t i=0;i<=amount;++i) {
		hashtab->pointer[i] = 0;
	}

	hashtab->seq_bkt = new uint32_t [len_genome-len_sed + 1];

	if ( NULL == hashtab->seq_bkt ) {
		fprintf(stderr, "Failed when applying for new space");
		exit(1);

	}

	uint32_t 	count = 1;
	// 读取临时哈希文件，排序并合并到最终的哈希索引文件中
	for(int i=0;i<NUM_FILE;i++) {
		//fprintf(stderr,"Read the %d file",i);
		seed *sed = new seed[count_line[i]];
		//fp[i]>>count_line[i]>>endl;
		fp[i].seekg(0,ios::beg);
		for(uint32_t j=0;j<count_line[i];j++) {
			//fp[i]>>sed[j].sed_value>>sed[j].sed_bkt_pos;
			fp[i].read((char*)(sed+j),sizeof(seed));

		}
		qsort(sed,count_line[i],sizeof(seed),comparator);
		//count_real_line = statsed(sed,count_line[i]);
		//fp[i].seekg(0,ios::beg);
		//fp[i]<<count_real_line<<endl;
		//fp[i].write((char *)&count_real_line,sizeof(uint32_t));

		num_pre = i << move_times;
		uint32_t temp_sed_value = sed[0].sed_value;
		uint32_t temp_sed_bkt_pos = sed[0].sed_bkt_pos;
		//fp[i] << sed[0].sed_value << " "<<sed[0].sed_bkt_pos << endl;
		
		//fp[i].write((char *)sed,sizeof(seed));
		hashtab->pointer[num_pre|temp_sed_value] = count;
		hashtab->seq_bkt[count - 1] = temp_sed_bkt_pos;
		++count;

		for (uint32_t j=1;j<count_line[i];++j) {
			if (sed[j].sed_value != temp_sed_value ) {
				//fp[i].write((char *)(sed+j),sizeof(seed));
				//fp[i] << sed[j].sed_value << " "<<sed[j].sed_bkt_pos << endl;
				hashtab->pointer[(num_pre|temp_sed_value)+1] = count;
				hashtab->pointer[num_pre|sed[j].sed_value] = count;
				hashtab->seq_bkt[count-1] = sed[j].sed_bkt_pos;
				temp_sed_value = sed[j].sed_value;
				temp_sed_bkt_pos = sed[j].sed_bkt_pos;
				++count;	
			} else {
				if (sed[j].sed_bkt_pos != temp_sed_bkt_pos){
					hashtab->seq_bkt[count-1] = sed[j].sed_bkt_pos;
					temp_sed_bkt_pos = sed[j].sed_bkt_pos;
					++count;
				}
			}								

		}

		hashtab->pointer[(num_pre|temp_sed_value)+1] = count;
		delete[] sed;
		fp[i].close();

	}
	//fprintf(stderr,"step 2");
	// remove tempory files and create hash file 删除临时文件，并创建最终的哈希文件
	for(int i=0;i<NUM_FILE;i++) {
		snprintf(num,4,"%d",i);
		strcpy(filepath,path);
		strcat(filepath,"/");
		strcat(filepath,num);
		remove(filepath);
	}

	FILE * fp_hash;
	strcpy(filepath,path);
	strcat(filepath,"/");
	strcat(filepath,"hash");
	
	if((fp_hash = fopen(filepath,"wb")) == NULL) {
		fprintf(stderr, "Failed when creating hash file, now exit");
		exit(1);
	} //else 
	//remove(filepath);
	//something may be wrong.
	fwrite(hashtab->pointer,sizeof(uint32_t),amount+1,fp_hash);
	fwrite(hashtab->seq_bkt,sizeof(uint32_t),count,fp_hash);
	fclose(fp_hash);

	if(NULL != hashtab->pointer)
		delete[] hashtab->pointer;

	/*
	for (uint32_t i=0; i<=len_genome - len_sed; ++i) {
	    if( NULL != hashtab->seq_bkt[i])
	        delete[] hashtab->seq_bkt[i]; 
	}
	*/
	if (NULL!=hashtab->seq_bkt)
	    delete[] hashtab->seq_bkt;
	


	if(NULL != hashtab)
		delete hashtab;
	
	return 1; // 返回1表示成功完成哈希文件的创建
}


