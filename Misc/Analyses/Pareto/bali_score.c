
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>

#define MAXNAMES 20
#define MAXSEQ 1000
#define MAXLEN 5000
#define MAXLINE 2048

#define Boolean short
#define TRUE 1
#define FALSE 0
#define EOS '\0'

#define INT_SCALE_FACTOR 100 /* Scaling factor to convert float to integer for profile scores */
#define NUMRES 32		/* max size of comparison matrix */

#define LEFT 1
#define RIGHT 2

#define NODE 0
#define LEAF 1


typedef struct node {		/* phylogenetic tree structure */
        struct node *left;
        struct node *right;
        struct node *parent;
        float dist;
        int  leaf;
        int order;
        char name[64];
} stree, *treeptr;

/* trees.c */
void guide_tree(FILE *tree,int nseqs);

int SeqGCGCheckSum(char *seq, int len);

Boolean keyword(char *line,char *code);
int getargs(char *line,char *args[],int max);
int parse_repart(char *str);
Boolean blankline(char *line);
void readmsf(FILE *fin);
void readref(FILE *fin);
void score_ref(void);
int countmsf(FILE *fin);
void columnscore(Boolean verbose);
void readref(FILE *fin);
int read_annotation(FILE *fin);
int countref(FILE *ifd);
void code_refseq(int i);
void code_seq(int i);
void aln_seq(int i);
void ref_gaps(int cutoff);

/* arrays for storing reference alignment */
char refnames[MAXSEQ][MAXNAMES+1];
char refseq_array[MAXSEQ][MAXLEN];
int refseq_code[MAXSEQ][MAXLEN];
int refseq_col[MAXLEN];
int refseqlength;
int refnseqs;

/* arrays for storing text alignment */
char names[MAXSEQ][MAXNAMES+1];
char seq_array[MAXSEQ][MAXLEN];
int seq_code[MAXSEQ][MAXLEN];
int seq_aln[MAXSEQ][MAXLEN];
int seqlength;
int seq_xref[MAXSEQ];
int nseqs;

Boolean verbose;

int refscore;
/* number of gaps in a column in the reference alignment */
int gap[MAXLEN];
/* cutoff for number of gaps allowed in a column in the reference alignment 
for column to be included in alignment score */
int cutoff;


/* core block data read in from annotation file */
int nblocks;
typedef struct Block {
	int s;
	int e;
	int nseqs;
} Block;

Block blocks[100];

/* alignment scores - mean square score, total column score and percentage
of correctly aligned pairs of residues */
float ms_score,tc_score;
float maxscore,maxpc_res;
float pc_res;

int main(int argc, char **argv)
{
	FILE *ifd,*tfd,*afd;
	int  err,i,j,ires,iseq=0;
	int ix;
	char t[MAXLINE+1];
	char seq[MAXLINE+1];
	char clen[MAXLINE+1];
	char cvar[MAXLINE+1];
	char cprog[MAXLINE+1];
	char method;
	Boolean eof,found;

	if(argc<3) {
		fprintf(stderr,"Usage: %s ref_aln test_aln [annot_file] [-v]\n",argv[0]);
		fprintf(stderr,"                where ref_aln       reference alignment in msf format \n");
		fprintf(stderr,"                      test_aln      test alignment in msf format \n");
		fprintf(stderr,"                      annot_file    optional BAliBASE annotation file\n");
		fprintf(stderr,"                      -v            verbose mode\n");
		return;
	}
/* open the reference aln file */

        if((ifd=fopen(argv[1],"r"))==NULL) {
            fprintf(stderr,"Cannot open reference aln file [%s]",argv[1]);
            return;
        }
/* open the test aln file */

        if((tfd=fopen(argv[2],"r"))==NULL) {
            fprintf(stderr,"Cannot open test aln file [%s]",argv[2]);
            return;
        }

/* open the annotation file */
/* plus a bit of trickery to decide if an annotation file is supplied */
/* set method="C" if no annotation file ie. we're all the columns in the
alignment (minus the columns with gaps in!) */

	verbose=FALSE;
	method='C';
	if(argc>=4)
	{
		if(argv[3][0]=='-')
			verbose=TRUE;
		else
		{
        		if((afd=fopen(argv[3],"r"))==NULL) {
            		fprintf(stderr,"Cannot open annotation file [%s]",argv[3]);
            		return;
			}
			method='B';
			if(argc>=5) verbose=TRUE;
		}
        }

/* read the reference alignment into refnames,refseq_array,refseqlength */
	refnseqs=countmsf(ifd);
	if(refnseqs==0) {
		fprintf(stderr,"Error: no sequences in %s\n",argv[1]);
		return;
	}
	readref(ifd);

/* read the test alignment into names, seq_array, seqlength */
	nseqs = countmsf(tfd);
	if(nseqs==0) {
		fprintf(stderr,"Error: no sequences in %s\n",argv[2]);
		return;
	}
	readmsf(tfd);

	if(nseqs != refnseqs) {
		fprintf(stderr,"Error: %d sequences in %s and %d in %s",refnseqs,argv[1],nseqs,argv[2]);
		return;
	}


/* if no annotation file supplied, only consider columns with less than 'cutoff' 
gaps. Here we use 20% of the number of sequences */
	if(method=='C')
	{
        	cutoff=(float)refnseqs*20.0/100.0;
        	if(cutoff<1) cutoff=1;
        	ref_gaps(cutoff);
	}

/* cross-reference the sequence refnames, in case they're not in the same order
in the reference and test alignment files */
	for(i=0;i<refnseqs;i++)
		seq_xref[i]=-1;
	for(i=0;i<refnseqs;i++)
	{
		found=FALSE;
		for(j=0;j<nseqs;j++)
		{
			if(strcasecmp(names[j],refnames[i])==0)
			{
				found=TRUE;
				seq_xref[j]=i;
				break;
			}
		}
		if(found==FALSE) {
			fprintf(stderr,"Error: sequence %s not found in test aln %s\n",refnames[i],argv[2]);
			return;
		}
	}

        fprintf(stdout,"\nComparing test alignment in %s\nwith reference alignment in %s\n",argv[2],argv[1]);
	if(method=='B')
	fprintf(stdout,"\nUsing core blocks defined in %s\n",argv[3]);

/* code the reference alignment - assign to each residue the number of the column it's in */
	for(i=0;i<refnseqs;i++)
		code_refseq(i);

/* read the annotation file */
	if(method=='B')
	{
		err=read_annotation(afd);
		if(err<0) fprintf(stderr,"No motif positions in %s: using all conserved columns\n",argv[4]);
	}

/* calculate the max score possible ie the score for the reference alignment */
	score_ref();

/* code the test alignment - look up each residue from the test alignment in the reference
alignment and assign the reference column number */
	for(i=0;i<nseqs;i++)
		if ( seq_xref[i] != -1) code_seq(i);

/* calculate the scores */
	columnscore(verbose);
	ms_score/=maxscore;
	pc_res/=maxpc_res;


        fprintf(stdout,"\n\tSP score= %.3f\n",pc_res);
        fprintf(stdout,"\n\tTC score= %.3f\n",tc_score);

	
}




void code_refseq(int seq)
{
	int i,j;

/* assign column no. in reference alignment to each residue */

	for(i=0;i<refseqlength;i++)
	{
               if(refseq_array[seq][i]=='-')
		{
			refseq_code[seq][i]=0;
		}
		else
		{
			refseq_code[seq][i]=i+1;
		}
	}
}




void ref_gaps(int cutoff)
{
	int i,j;

/* find columns with gaps in the reference sequnce - set refseq_col[]=0
if gaps, =nseqs otherwise */

	for(i=0;i<refseqlength;i++)
	{
                gap[i]=0;
                for(j=0;j<refnseqs;j++)
                        if(refseq_array[j][i]=='-')
                        {
                                gap[i]++;
                        }
	}

	for(i=0;i<refseqlength;i++)
	{
		if(gap[i]>=cutoff) refseq_col[i]=0;
		else refseq_col[i]=refnseqs;
	}
}

void score_ref(void)
{
	int i,j;

/* calculates the maximum score possible for an alignment */
	
	for(i=0;i<refseqlength;i++)
	{
		if(refseq_col[i]>1) 
		{
			maxscore+=refseq_col[i]*refseq_col[i];
			maxpc_res+=refseq_col[i]*(refseq_col[i]-1)/2.0;
		}
	}
}

void code_seq(int seq)
{
	int i,j,ix;

/* find the first residue in the reference sequence */
	ix=0;
	for(j=0;j<refseqlength;j++)
		if(refseq_array[seq_xref[seq]][j]!='-')
		{
			ix=refseq_code[seq_xref[seq]][j];
			break;
		}

	for(i=0;i<seqlength;i++)
	{
		if(seq_array[seq][i]=='-')
		{
			seq_code[seq_xref[seq]][i]=0;
		}
		else
		{
			if (refseq_col[ix-1]>seq_xref[seq]) seq_code[seq_xref[seq]][i]=ix;
			for(j+=1;j<refseqlength;j++)
				if(refseq_array[seq_xref[seq]][j]!='-')
				{
					ix=refseq_code[seq_xref[seq]][j];
					break;
				}
		}
	}
}

void columnscore(Boolean verbose)
{
	int i,j,k;
	int iseq,ncols;
	int scores[MAXSEQ];
	int index[MAXSEQ];
	int n;
	float colscore[MAXLEN];
	float colscore1[MAXLEN];
	int pc;
	int nblocks;
	int blen=30;
	Boolean found;

	ms_score=tc_score=ncols=0;
	for(j=0;j<seqlength;j++)
	{
		colscore[j]=0;
		colscore1[j]=0;
	}
	for(i=0;i<seqlength;i++)
	{
		for(j=0;j<nseqs;j++)
			scores[j]=0;
		n=0;
		for(j=0;j<nseqs;j++)
		{
			if(seq_code[j][i]!=0)
			{
				found=FALSE;
				for(k=0;k<n;k++)
				{
					if(seq_code[j][i]==index[k])
					{
						scores[k]++;
						found=TRUE;
						break;
					}
				}
				if(found==FALSE)
				{
					scores[n]=1;
					index[n]=seq_code[j][i];
					n++;
				}
			}
		}
		for(j=0;j<nseqs;j++)
		{
			if(scores[j]>1)
			{
				colscore[i]+=scores[j]*scores[j];
				pc_res+=scores[j]*(scores[j]-1)/2.0;
			}
		}
/* count 1 for each column */
		for(j=0;j<nseqs;j++)
			if(seq_code[0][i]>0 && scores[j]>=refseq_col[seq_code[0][i]-1])
				colscore1[i]=1;
		if (seq_code[0][i]!=0) ncols++;

		ms_score+=colscore[i];
		tc_score+=colscore1[i];
	}

	if(verbose)
	{
		nblocks=seqlength/blen;
                fprintf(stdout,"\n\n");
		for(k=0;k<=nblocks;k++) {
        		for(i=0;i<nseqs;i++) {
                		fprintf(stdout,"%10s",names[i]);
               			for(j=0;j<blen && k*blen+j<seqlength;j++) {
                       			fprintf(stdout,"   %c",seq_array[i][k*blen+j]);
                		}
                		fprintf(stdout,"\n");
			}
                	fprintf(stdout,"\n          ");
               		for(j=0;j<blen && k*blen+j<seqlength;j++) {
                       		fprintf(stdout,"% 4.0f",colscore[k*blen+j]);
                	}
                	fprintf(stdout,"\n\n");
        	}

	}

/* count 1 for each column */
	tc_score/=(float)ncols;
}


void code_blockseq(int seq)
{
	int i,j;

	for(i=0,j=1;i<seqlength;i++)
	{
		if(seq_array[seq][i]=='-')
		{
			seq_code[seq][i]=0;
		}
		else
		{
			seq_code[seq][i]=j++;
		}
	}
}



int read_annotation(FILE *fin)
{
	int i,j,k,s,e;
	char c,line[MAXLINE+1];
	char t[MAXLINE+1];
	int  numargs;
	char *args[101];
	Boolean g;

	nblocks=0;
        while(fgets(line,MAXLINE+1,fin)) {
		if(keyword(line,"BPOS")) {
			nblocks=0;
			numargs = getargs(line,args,101);
			for(i=1;i<numargs;i++)
			{
				s=sscanf(args[i],"%d%c%d",&blocks[nblocks].s,&c,&blocks[nblocks].e);
				if(s==1) {
					blocks[nblocks].e=blocks[nblocks].s;
				}
				else if(s!=3) {
					break;
				}
				blocks[nblocks].nseqs=refnseqs;
				nblocks++;
			}
			break;
		}
	}
	
        while(fgets(line,MAXLINE+1,fin)) {
		if(keyword(line,"BREPART")) {
			numargs = getargs(line,args,101);
			for(i=1;i<numargs;i++)
			{
				blocks[i-1].nseqs = 0;
				sscanf(args[i],"%s",t);
				blocks[i-1].nseqs = parse_repart(t);
			}
		}
	}

	for(i=0;i<refseqlength;i++)
		refseq_col[i]=0;

	for(i=0;i<nblocks;i++) {
		if(blocks[i].nseqs > 0) {
			for(j=0, k=0;j<refseqlength;j++)
			{
				if(isalpha(refseq_array[0][j])) 
					k++;
				if(k==blocks[i].s) {
				s=j;
				break;
				}
			}
			for(;j<refseqlength;j++)
			{
				if(isalpha(refseq_array[0][j])) 
					k++;
				if(k==blocks[i].e) {
				e=j+1;
				break;
				}
			}
			for(j=s;j<=e;j++)
				refseq_col[j]=blocks[i].nseqs;
		}
	}

}

Boolean keyword(char *line,char *code)
{
        int i;
        char key[MAXLINE];

        for(i=0;!isspace(line[i]) && line[i]!=EOS;i++)
                key[i]=line[i];
        key[i]=EOS;
        return( strcmp(key,code) == 0 );
}

int getargs(char *line,char *args[],int max)
{

        char    *inptr;
        int     i;

        inptr=line;
        for (i=0;i<=max;i++)
        {
                if ((args[i]=strtok(inptr," \t\n"))==NULL)
                        break;
                inptr=NULL;
        }
        if (i>max)
        {
                fprintf(stdout,"Too many args\n");
                return(0);
        }

        return(i);
}

int parse_repart(char *str)
{
	int len,s,e,i,ix,tot;
	char t[MAXLINE+1];

/* string contains a number of substrings of format n<Hn or n, separated by + */
	tot=0;
	len=strlen(str);
	for(ix=0;ix<len;)
	{
		for(;ix<len && !isalnum(str[ix]);ix++);
		s=ix;
		for(;ix<len && str[ix]!='+';ix++);
		e=ix;
		for(i=s;i<e;i++)
			t[i-s]=str[i];
		t[i-s]='\0';
/* now we have a string in t[] with format n<Hn or n */
		if(toupper(t[0])=='H') {
			sscanf(&t[1],"%d",&i);
			if(i>0 && i<=refnseqs) tot += i;
			else fprintf(stderr,"WARNING: bad format in repartition %s\n",str);
		}
		else {
			sscanf(t,"%d",&i);
			if(i>0 && i<=refnseqs) tot += i;
			else fprintf(stderr,"WARNING: bad format in repartition %s\n",str);
		}
	}

	return tot;
}

Boolean blankline(char *line)
{
	int i;

        for(i=0;line[i]!='\n' && line[i]!=EOS;i++) {
                if( isdigit(line[i]) ||
                    isspace(line[i]) ||
                    (line[i] == '*') ||
                    (line[i] == ':') ||
                    (line[i] == '.'))
                        ;
                else
                        return FALSE;
        }
        return TRUE;
}

void readref(FILE *fin)
{
        static char line[MAXLINE+1];
        char seq[MAXLEN+1];
        int len,seqno,i,j,k;
        unsigned char c;

for(seqno=0;seqno<refnseqs;seqno++)
{

        fseek(fin,0,0);                 /* start at the beginning */

        len=0;                         /* initialise length to zero */
        for(i=0;;i++) {
                if(fgets(line,MAXLINE+1,fin)==NULL) return; /* read the title*/
                if(line[0]=='/' && line[1]=='/') break;
 
        }

        while (fgets(line,MAXLINE+1,fin) != NULL) {
                if(!blankline(line)) {

                        for(i=0;i<seqno;i++) fgets(line,MAXLINE+1,fin);
                        for(j=0;j<=strlen(line);j++) if(line[j] != ' ') break;
                        for(k=j;k<=strlen(line);k++) if(line[k] == ' ') break;
                        strncpy(refnames[seqno],line+j,k-j);
                        refnames[seqno][k-j]='\0';
                        refnames[seqno][MAXNAMES]='\0';

                        for(i=k;i<=MAXLINE;i++) {
                                c=line[i];
                                if(c == '.' || c == '~' ) c = '-';
                                if(c == '*') c = 'X';
                                if(c == '\n' || c == EOS) break; /* EOL */
                                if(isalpha(c) || c=='-') seq[len++]=c;
                        }

                        for(i=0;;i++) {
                                if(fgets(line,MAXLINE+1,fin)==NULL) break;
                                if(blankline(line)) break;
                        }
                }
        }
	strcpy(refseq_array[seqno],seq);
	if(seqno==0) refseqlength=len;
}
}

void readmsf(FILE *fin)
{
        static char line[MAXLINE+1];
        char seq[MAXLEN+1];
        int len,seqno,i,j,k;
        unsigned char c;

for(seqno=0;seqno<nseqs;seqno++)
{

        fseek(fin,0,0);                 /* start at the beginning */

        len=0;                         /* initialise length to zero */
        for(i=0;;i++) {
                if(fgets(line,MAXLINE+1,fin)==NULL) return; /* read the title*/
                if(line[0]=='/' && line[1]=='/') break;
 
        }

        while (fgets(line,MAXLINE+1,fin) != NULL) {
                if(!blankline(line)) {

                        for(i=0;i<seqno;i++) fgets(line,MAXLINE+1,fin);
                        for(j=0;j<=strlen(line);j++) if(line[j] != ' ') break;
                        for(k=j;k<=strlen(line);k++) if(line[k] == ' ') break;
                        strncpy(names[seqno],line+j,k-j);
                        names[seqno][k-j]='\0';
                        names[seqno][MAXNAMES]='\0';

                        for(i=k;i<=MAXLINE;i++) {
                                c=line[i];
                                if(c == '.' || c == '~' ) c = '-';
                                if(c == '*') c = 'X';
                                if(c == '\n' || c == EOS) break; /* EOL */
                                if(isalpha(c) || c=='-') seq[len++]=c;
                        }

                        for(i=0;;i++) {
                                if(fgets(line,MAXLINE+1,fin)==NULL) break;
                                if(blankline(line)) break;
                        }
                }
        }
	strcpy(seq_array[seqno],seq);
	if(seqno==0) seqlength=len;
}
}

int countmsf(FILE *fin)
{
/* count the number of sequences in a PILEUP alignment file */

        char line[MAXLINE+1];
        int  lnseqs;

        while (fgets(line,MAXLINE+1,fin) != NULL) {
                if(line[0]=='/' && line[1]=='/') break;
        }

        while (fgets(line,MAXLINE+1,fin) != NULL) {
                if(!blankline(line)) break;             /* Look for next non- */
        }                                               /* blank line */
        lnseqs = 1;

        while (fgets(line,MAXLINE+1,fin) != NULL) {
                if(blankline(line)) return lnseqs;
                lnseqs++;
        }

        return 0; /* if you got to here-funny format/no seqs.*/
}
