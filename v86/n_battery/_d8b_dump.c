
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define SFA_IMPLEMENTATION
#include "../../sfa/format/sfa.h"
int main(int c,char**v){SFA*s=sfa_open(v[1]);if(!s)return 1;
 size_t nb=0;for(uint32_t i=0;i<s->n_columns;i++)nb+=s->N_total*sfa_dtype_size[s->columns[i].dtype];
 void*b=malloc(nb); int fr=atoi(v[2]); if(sfa_read_frame(s,fr,b)!=0)return 2;
 fprintf(stderr,"%u %.10g %u\n",s->Nx,s->Lx,s->total_frames);
 const char*want[12]={"phi_x","phi_y","phi_z","phiim_x","phiim_y","phiim_z",
   "phi_vx","phi_vy","phi_vz","phiim_vx","phiim_vy","phiim_vz"};
 for(int w=0;w<12;w++){uint64_t off=0;int found=-1;
  for(uint32_t i=0;i<s->n_columns;i++){if(!strncmp(s->columns[i].name,want[w],32)){found=i;break;}
   off+=s->N_total*sfa_dtype_size[s->columns[i].dtype];}
  if(found<0){fprintf(stderr,"missing %s\n",want[w]);return 3;}
  fwrite((char*)b+off,4,s->N_total,stdout);} return 0;}
