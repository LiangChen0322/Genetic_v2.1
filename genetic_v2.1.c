/******************************************************************/
/* 基于基本遗传算法的函数最优化 SGA.C */
/* A Function Optimizer using Simple Genetic Algorithm */
/* developed from the Pascal SGA code presented by David E.Goldberg */
//******************************************************************/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "stdlib.h"
/* 全局变量 */
struct individual /* 个体*/
{
	unsigned *chrom; /* 染色体 */
	double fitness; /* 个体适应度*/
	double varible; /* 个体对应的变量值*/ 
	int xsite; /* 交叉位置 */
	int parent[2]; /* 父个体 */
	/* 特定数据指针变量 */
};
struct bestever /* 最佳个体*/
{
	unsigned *chrom; /* 最佳个体染色体*/
	double fitness; /* 最佳个体适应度 */
	double varible; /* 最佳个体对应的变量值 */
	int generation; /* 最佳个体生成代 */
};
struct individual *oldpop; /* 当前代种群 */
struct individual *newpop; /* 新一代种群 */
struct bestever bestfit; /* 最佳个体 */
double sumfitness; /* 种群中个体适应度累计 */
double max; /* 种群中个体最大适应度 */
double avg; /* 种群中个体平均适应度 */
double min; /* 种群中个体最小适应度 */
float pcross; /* 交叉概率 */
float pmutation; /* 变异概率 */
int popsize; /* 种群大小 */
int lchrom; /* 染色体长度*/
int chromsize; /* 存储一染色体所需字节数 */
int gen; /* 当前世代数 */
int maxgen; /* 最大世代数 */
int run; /* 当前运行次数 */
int maxruns; /* 总运行次数 */
int printstrings; /* 输出染色体编码的判断，0 -- 不输出, 1 -- 输出 */
int nmutation; /* 当前代变异发生次数 */
int ncross; /* 当前代交叉发生次数 */

/* 随机数发生器使用的静态变量 */
static double oldrand[55];
static int jrand;
static double rndx2;
static int rndcalcflag;
/* 输出文件指针 */
FILE *outfp ;
/* 函数定义 */
void advance_random();
int flip(float);
int rnd(int, int);
void randomize();
float randomperc();
void warmup_random(float);
void initialize(),initdata(),initpop();
void initreport(),generations(),initmalloc();
void freeall(),nomemory(char *string),report();
void writepop(),writechrom(unsigned *chrom);
void preselect_p();
void statistics(struct individual *pop);
void title(),repchar (FILE *,char *,int);

int select_p();
void objfunc(struct individual *critter);
int crossover (unsigned *, unsigned *, unsigned *, unsigned *);
void mutation(unsigned *);

void initialize() /* 遗传算法初始化 */
{
	/* 键盘输入遗传算法参数 */
	initdata();
	/* 确定染色体的字节长度 */
	chromsize = (lchrom/(8*sizeof(unsigned)));
	if(lchrom%(8*sizeof(unsigned)))  
		chromsize++;
	/*分配给全局数据结构空间 */
	initmalloc();
	/* 初始化随机数发生器 */
	randomize();
	/* 初始化全局计数变量和一些数值*/
	nmutation = 0;
	ncross = 0;
	bestfit.fitness = 0.0;
	bestfit.generation = 0;
	/* 初始化种群，并统计计算结果 */
	initpop();
	statistics(oldpop);
	initreport();
}

void initdata() /* 遗传算法参数输入 */
{
	char answer[10];
	printf("种群大小为(20-100)：");
	scanf("%d", &popsize);
	if((popsize%2) != 0)
	{ 
		printf("种群大小已设置为偶数\n");
		popsize++;
	};

	printf("染色体长度(8-40)：");
	scanf("%d", &lchrom);

	printf("是否输出染色体编码?(yes or no)：");
	printstrings=1;
	scanf("%s", answer);
	if(strncmp(answer,"n",1) == 0)
		printstrings = 0;

	printf("最大世代数(100-300)：");
	scanf("%d", &maxgen);

	printf("交叉率(0.2-0.9)：");
	scanf("%f", &pcross);

	printf("变异率(0.01-0.1)：");
	scanf("%f", &pmutation);
}

/* 随机初始化种群 */
void initpop() 
{
	int j, j1, k, stop;
	unsigned mask = 1;
	for(j = 0; j < popsize; j++) {
		for(k = 0; k < chromsize; k++) {
			oldpop[j].chrom[k] = 0;
			if(k == (chromsize-1)) {
				stop = lchrom - (k*(8*sizeof(unsigned)));
			} else {
				stop =8*sizeof(unsigned);
			}
			for(j1 = 1; j1 <= stop; j1++) {
				oldpop[j].chrom[k] = oldpop[j].chrom[k]<<1;
				if(flip(0.5)) {
					oldpop[j].chrom[k] = oldpop[j].chrom[k]|mask;
				}
			}
		}
		oldpop[j].parent[0] = 0; /* 初始父个体信息 */
		oldpop[j].parent[1] = 0;
		oldpop[j].xsite = 0;
		objfunc(&(oldpop[j])); /* 计算初始适应度*/
	}
}

void initreport() /* 初始参数输出 */
{
	fprintf(outfp," 基本遗传算法参数\n");
	fprintf(outfp," -------------------------------------------------\n");
	fprintf(outfp," 种群大小(popsize) = %d\n",popsize);
	fprintf(outfp," 染色体长度(lchrom) = %d\n",lchrom);
	fprintf(outfp," 最大进化代数(maxgen) = %d\n",maxgen);
	fprintf(outfp," 交叉概率(pcross) = %f\n", pcross);
	fprintf(outfp," 变异概率(pmutation) = %f\n", pmutation);
	fprintf(outfp," -------------------------------------------------\n");
	fflush(outfp);
}

void generations()
{
	int mate1, mate2, jcross, j = 0;
	/* 每代运算前进行预选 */
	preselect_p();
	/* 选择, 交叉, 变异 */
	do
	{
		/* 挑选交叉配对 */
		mate1 = select_p();
		mate2 = select_p();
		/* 交叉和变异 */
		jcross = crossover(oldpop[mate1].chrom, oldpop[mate2].chrom, newpop[j].chrom, newpop[j+1].chrom);
		mutation(newpop[j].chrom);
		mutation(newpop[j+1].chrom);
		/* 解码, 计算适应度 */
		objfunc(&(newpop[j]));
		/*记录亲子关系和交叉位置 */
		newpop[j].parent[0] = mate1+1;
		newpop[j].xsite = jcross;
		newpop[j].parent[1] = mate2+1;
		objfunc(&(newpop[j+1]));
		newpop[j+1].parent[0] = mate1+1;
		newpop[j+1].xsite = jcross;
		newpop[j+1].parent[1] = mate2+1;
		j = j + 2;
	}
	while(j < (popsize-1));

}

void initmalloc() /*为全局数据变量分配空间 */
{
	unsigned nbytes;

	int j;
	/* 分配给当前代和新一代种群内存空间 */
	nbytes = popsize*sizeof(struct individual);
	if((oldpop = (struct individual *)malloc(nbytes)) == NULL)
		nomemory("oldpop");
	if((newpop = (struct individual *)malloc(nbytes)) == NULL)
		nomemory("newpop");
	/* 分配给染色体内存空间 */
	nbytes = chromsize*sizeof(unsigned);
	for(j = 0; j < popsize; j++) {
		if((oldpop[j].chrom = (unsigned *)malloc(nbytes)) == NULL)
			nomemory("oldpop chromosomes");
		if((newpop[j].chrom = (unsigned *)malloc(nbytes)) == NULL)
			nomemory("newpop chromosomes");
	}
	if((bestfit.chrom = (unsigned *)malloc(nbytes)) == NULL)
		nomemory("bestfit chromosome");

}

void freeall() /* 释放内存空间 */
{
	int i;
	for(i = 0; i < popsize; i++) {
		free(oldpop[i].chrom);
		free(newpop[i].chrom);
	}
	free(oldpop);
	free(newpop);
	free(bestfit.chrom);
}

void nomemory(char *string) /* 内存不足，退出*/
{
	printf("malloc: out of memory making %s!!\n",string);
	exit(-1);
}

void report() /* 输出种群统计结果 */
{

	extern void writepop();

	if(printstrings == 1) {
        fprintf(outfp,"模拟计算统计报告 \n");
		fprintf(outfp, "世代数 %3d\n", gen);
        fprintf(outfp,"个体  染色体编码  适应度");
        fprintf(outfp," 父个体 交叉位置 染色体编码 适应度\n");
        writepop();
    }
	fprintf(outfp,"第 %d 代统计: \n",gen);
	fprintf(outfp,"总交叉操作次数 = %d, 总变异操作数 = %d\n",ncross,nmutation);
	fprintf(outfp," 最小适应度：%f 最大适应度：%f 平均适应度 %f\n", min,max,avg);
	fprintf(outfp," 迄今发现最佳个体 => 所在代数： %d \n", bestfit.generation);
	fprintf(outfp," 最佳个体适应度：%f \n", bestfit.fitness);
	fprintf(outfp," 最佳个体染色体：\n");
	writechrom((&bestfit)->chrom);
	fprintf(outfp,"最佳个体对应的变量值: %f\n", bestfit.varible);
	fprintf(outfp," -------------------------------------------------\n");

}

void writepop()
{
	struct individual *pind;
	int j;

	for(j=0; j<popsize; j++) {
		fprintf(outfp,"%3d" ,j+1);
		/* 当前代个体 */
		pind = &(oldpop[j]);
		writechrom(pind->chrom);
		fprintf(outfp,"%8f", pind->fitness);
		/* 新一代个体 */
		pind = &(newpop[j]);
		fprintf(outfp,"(%2d,%2d) %2d ",
			pind->parent[0], pind->parent[1], pind->xsite);
		writechrom(pind->chrom);
		fprintf(outfp," %8f\n", pind->fitness);
	}
}

/* 输出染色体编码 */
void writechrom(unsigned *chrom) 
{
	int j, k, stop;
	unsigned mask = 1, tmp;

	for (k = 0; k < chromsize; k++) {
		tmp = chrom[k];
		if (k == (chromsize-1))
			stop = lchrom - (k*(8*sizeof(unsigned)));
		else
			stop =8*sizeof(unsigned);
		for (j = 0; j < stop; j++) {
			if(tmp&mask)
				fprintf(outfp,"1");
			else
				fprintf(outfp,"0");
			tmp = tmp>>1;
		}
	}
}

void preselect_p()
{
	int j;
	sumfitness = 0;
	for (j = 0; j < popsize; j++) 
		sumfitness += oldpop[j].fitness;
}

int select_p() /*轮盘赌选择*/
{
	extern float randomperc();
	float sum, pick;
	int i;
	
	pick = randomperc();
	sum = 0;
	if (sumfitness != 0) {
		for(i = 0; (sum < pick) && (i < popsize); i++)
			sum += oldpop[i].fitness/sumfitness;
	}
	else
		i = rnd(1,popsize);
	return(i-1);
}

void statistics(struct individual *pop) /* 计算种群统计数据 */
{
	int i, j;
	sumfitness = 0.0;
	min = pop[0].fitness;
	max = pop[0].fitness;
	/* 计算最大、最小和累计适应度 */
	for (j = 0; j < popsize; j++) {
		sumfitness = sumfitness + pop[j].fitness; 
		if(pop[j].fitness > max) 
			max = pop[j].fitness; 
		if(pop[j].fitness < min) 
			min = pop[j].fitness; 
		/* new global best-fit individual */
		if (pop[j].fitness > bestfit.fitness) {
			for (i = 0; i < chromsize; i++)
				bestfit.chrom[i] = pop[j].chrom[i];
			bestfit.fitness = pop[j].fitness;
			bestfit.varible = pop[j].varible; 
			bestfit.generation = gen;
		}
	}
	/* 计算平均适应度 */
	avg = sumfitness/popsize;
}

/* 计算适应度函数值 */
void objfunc(struct individual *critter)
{
	unsigned mask = 1;
	unsigned bitpos;
	unsigned tp;
	double bitpow ;
	int j, k, stop;
	
	critter->varible = 0.0;
	for (k = 0; k < chromsize; k++) {
		if (k == (chromsize-1))
			stop = lchrom-(k*(8*sizeof(unsigned)));
		else
			stop =8*sizeof(unsigned);
		tp = critter->chrom[k];
		for (j = 0; j < stop; j++) {
			bitpos = j + (8*sizeof(unsigned))*k;
			if((tp&mask) == 1) {
				bitpow = pow(2.0,(double) bitpos);
				critter->varible = critter->varible + bitpow;
			}
			tp = tp>>1;
		}
	}
	critter->varible = -1+critter->varible*3/(pow(2.0,(double)lchrom)-1);
	//critter->fitness = critter->varible*sin(critter->varible*10*atan(1.0)*4)+2.0;
	critter->fitness = critter->varible*2-(critter->varible*critter->varible)+2.0;
}

/*变异操作*/
void mutation(unsigned *child) 
{
	int j, k, stop;
	unsigned mask, temp = 1;
	for(k = 0; k < chromsize; k++)
	{
		mask = 0;
		if(k == (chromsize-1))
			stop = lchrom - (k*(8*sizeof(unsigned)));
		else
			stop = 8*sizeof(unsigned);
		for(j = 0; j < stop; j++) {
			if(flip(pmutation)) {
				mask = mask|(temp<<j);
				nmutation++;
			}
		}
		child[k] = child[k]^mask;
	}
}
/* 由两个父个体交叉产生两个子个体 */
int crossover (unsigned *parent1, unsigned *parent2, 
				unsigned *child1, unsigned *child2)
{
	int j, jcross, k;
	unsigned mask, temp;

	if (flip(pcross)) {
		jcross = rnd(1 ,(lchrom - 1));	/* Cross between 1 and l-1 */
		ncross++;
		for(k = 1; k <= chromsize; k++) {
			if(jcross >= (k*(8*sizeof(unsigned)))) {
				child1[k-1] = parent1[k-1];
				child2[k-1] = parent2[k-1];
			} else if ((jcross < (k*(8*sizeof(unsigned)))) 
						&& (jcross > ((k-1)*(8*sizeof(unsigned))))) {
				mask = 1;
				for(j = 1; j <= (jcross-1-((k-1)*(8*sizeof(unsigned)))); j++) {
					temp = 1;
					mask = mask<<1;
					mask = mask|temp;
				}
				child1[k-1] = (parent1[k-1]&mask)|(parent2[k-1]&(~mask));
				child2[k-1] = (parent1[k-1]&(~mask))|(parent2[k-1]&mask);
			} else {
				child1[k-1] = parent2[k-1];
				child2[k-1] = parent1[k-1];
			}
		}
	} else {
		for (k = 0; k < chromsize; k++) {
			child1[k] = parent1[k];
			child2[k] = parent2[k];
		}
		jcross = 0;
	}
	return(jcross);
}

/* 产生55个随机数 */
void advance_random() 
{
	int j1;
	double new_random;
	
	for(j1 = 0; j1 < 24; j1++) {
		new_random = oldrand[j1] - oldrand[j1+31];
		if(new_random < 0.0) 
			new_random = new_random + 1.0;
		oldrand[j1] = new_random;
	}
	for(j1 = 24; j1 < 55; j1++) {
		new_random = oldrand [j1] - oldrand [j1-24];
		if(new_random < 0.0) 
			new_random = new_random + 1.0;
		oldrand[j1] = new_random;
	}
}

int flip(float prob) /* 以一定概率产生0或1 */
{
	float randomperc();
	if(randomperc() <= prob)
		return(1);
	else
		return(0);
}

/* 设定随机数种子并初始化随机数发生器 */
void randomize()
{
	float randomseed;
	int j1;
	for(j1=0; j1<=54; j1++)
		oldrand[j1] = 0.0;
	jrand=0;
	do
	{

		printf("随机数种子[0-1]:");
		scanf("%f", &randomseed);
	}
	while((randomseed < 0.0) || (randomseed > 1.0));
	warmup_random(randomseed);
}

/*与库函数random()作用相同, 产生[0,1]之间一个随机数 */
float randomperc()
{
	jrand++;
	if(jrand >= 55) {
		jrand = 1;
		advance_random();
	}
	return((float) oldrand[jrand]);
}

int rnd(int low, int high) /*在整数low和high之间产生一个随机整数*/
{
	int i;
	float randomperc();
	if(low >= high)
		i = low;
	else
	{
		i = (randomperc() * (high - low + 1)) + low;
		if(i > high) i = high;
	}
	return(i);
}

void warmup_random(float random_seed) /* 初始化随机数发生器*/
{
	int j1, ii;
	double new_random, prev_random;
	new_random = 0.000000001;
	prev_random = random_seed;
	for(j1 = 1 ; j1 <= 54; j1++)
	{
		ii = (21*j1)%54;
		oldrand[ii] = new_random;
		new_random = prev_random-new_random;
		if(new_random<0.0) 
			new_random = new_random + 1.0;
		prev_random = oldrand[ii];
	}
	advance_random();
	advance_random();
	advance_random();
	jrand = 0;
}

 main() /* 主程序 */
{
	struct individual *temp;
	outfp=fopen("output8.txt","w");
	printf("输入遗传算法执行次数(1-5):");
	scanf("%d",&maxruns);
	for(run=1; run<=maxruns; run++)
	{
		initialize();
		for(gen=0; gen<maxgen; gen++)
		{
			fprintf(outfp,"\n第 %d /%d 次运行: 当前代为 %d, 共 %d 代\n", run,maxruns,gen,maxgen);
			/* 产生新一代 */
			generations();
			/* 计算新一代种群的适应度统计数据 */
			statistics(newpop);
			/* 输出新一代统计数据 */
			report();
			temp = oldpop;
			oldpop = newpop;
			newpop = temp;
		}
		freeall();
	}
}