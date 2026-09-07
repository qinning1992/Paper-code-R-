#!/usr/bin/env Rscript
# Complete S13-S15, reconstructing the models described in the source supplement.
# Usage: Rscript Complete_S13_S15_HS_sensitivity.R <project_dir> <output_dir> [patient_xlsx]
# Raw input is read-only. Pathology_Group==1 defines HS (39/32), NOT the 40/31 HS flag.
# Scope: 19 core outcomes only (auxiliary WM SS removed per revision).
# Original S13-S15 executable/complete exports were unavailable. Effect sizes and
# nested F tests can be compared exactly up to published rounding; old Monte Carlo
# P/q values cannot be exactly reproduced without the old RNG stream.
suppressPackageStartupMessages({library(readxl);library(dplyr)})
invisible(Sys.setlocale('LC_CTYPE','English_United States.utf8'))
args=commandArgs(TRUE);stopifnot(length(args)>=2)
project=args[1];out=args[2];input=if(length(args)>2)args[3] else file.path(project,'code and data/df_raw_clean.xlsx')
rd=file.path(out,'02_results');vd=file.path(out,'06_validation')
dir.create(rd,recursive=TRUE,showWarnings=FALSE);dir.create(vd,recursive=TRUE,showWarnings=FALSE)
B=5000L;SEED=2026L
d=as.data.frame(read_excel(input));d=d[!is.na(d$side),]
stopifnot(nrow(d)==71,all(d$side %in% c(0,1)))
d$sex=factor(d$sex);d$HS_status=factor(ifelse(d$Pathology_Group==1,'HS','non-HS'))
stopifnot(sum(d$HS_status=='HS')==39,sum(d$HS_status=='non-HS')==32)
for(c in c('Duration_Dis','Edu_year','eTIV','Age'))d[[c]]=as.numeric(d[[c]])
original=read.csv(file.path(project,'analysis_outputs/new_sensitivity_analyses/TableS8_HS_status_stratified_sensitivity.csv'),check.names=FALSE,fileEncoding='UTF-8')
sc=unique(original[,c('Scale','Scale_Display','Is_Normed')]);sc=sc[sc$Scale!='加工记忆量表分',,drop=FALSE];stopifnot(nrow(sc)==19)
sc$Scope='Core battery'  # WM SS removed; 19 core outcomes only
rois=unique(original[,c('Family','ROI','ROI_Display')]);stopifnot(nrow(rois)==20)
defs=list(AnteriorDorsalGroup=c('AV','LD'),MediodorsalNuclei=c('MDm','MDl'),IntralaminarNuclei=c('CL','CM','CeM','Pf','Pc'),PulvinarLPComplex=c('Pu_Total','L_Sg','LP'),VentralNuclei=c('VA','VAmc','VLa','VLp','VPL','VM'),MidlineNuclei=c('MV_Re','Pt'),GeniculateNuclei=c('LGN','MGN'))
segments=c('Whole_hippocampal_head','Whole_hippocampal_body','Hippocampal_tail')
for(s in segments)for(o in c('ipsi','contra')){
 left=if(o=='ipsi')d$side==1 else d$side==0
 d[[paste0(o,'_hp_',s)]]=ifelse(left,d[[paste0('lh_',s)]],d[[paste0('rh_',s)]])
}
for(g in names(defs))for(o in c('ipsi','contra')){
 left=if(o=='ipsi')d$side==1 else d$side==0
 d[[paste0(o,'_thalgrp_',g)]]=ifelse(left,rowSums(d[,paste0('Left_',defs[[g]])]),rowSums(d[,paste0('Right_',defs[[g]])]))
}
d$ipsi_whole_hippocampus=ifelse(d$side==1,d[['Left-Subcort-Hippocampus']],d[['Right-Subcort-Hippocampus']])
stopifnot(max(abs(d$ipsi_whole_hippocampus-rowSums(d[,paste0('ipsi_hp_',segments)])))<1e-5)
# Include the overall ipsilateral hippocampus required by the reviewer's question.
# Treat all seven hippocampal measures as one correction family per cognitive outcome.
rois=rbind(data.frame(Family='Hippocampal measures',ROI='ipsi_whole_hippocampus',ROI_Display='Overall ipsilateral hippocampus'),rois)
rois$Family[rois$Family=='Hippocampal segments']='Hippocampal measures'
stopifnot(nrow(rois)==21)
covars=function(normed,hs=FALSE){z=c('sex','Duration_Dis','Edu_year','eTIV');if(!normed)z=c(z,'Age');if(hs)z=c(z,'HS_status');z}
seed_for=function(key){v=utf8ToInt(enc2utf8(key));as.integer((SEED+sum(as.double(v)*seq_along(v)))%%2147483646)+1L}
corr=function(data,x,y,cv,key){
 a=data[complete.cases(data[,c(x,y,cv)]),];X=model.matrix(reformulate(cv),a);Q=qr(X)
 stopifnot(Q$rank==ncol(X),nrow(a)>=15)
 rx=qr.resid(Q,as.numeric(a[[x]]));ry=qr.resid(Q,as.numeric(a[[y]]));r=cor(rx,ry)
 n=nrow(a);k=Q$rank-1L;df=n-k-2L
 set.seed(seed_for(key));idx=replicate(B,sample.int(n));xp=as.numeric(scale(rx));yp=as.numeric(scale(ry))
 rp=as.numeric(crossprod(xp,matrix(yp[idx],nrow=n)))/(n-1)
 p=(1+sum(abs(rp)>=abs(r)-1e-12))/(B+1)
 ci=tanh(atanh(r)+c(-1,1)*qnorm(.975)/sqrt(n-k-3))
 c(n=n,Covariate_df=k,Residual_df=df,r=r,CI_low=ci[1],CI_high=ci[2],P_perm=p)
}
S13=list();S15=list();n=0
for(i in seq_len(nrow(rois))){
 for(j in seq_len(nrow(sc))){
  n=n+1;x=rois$ROI[i];y=sc$Scale[j]
  meta=data.frame(Family=rois$Family[i],ROI=x,ROI_Display=rois$ROI_Display[i],Scale=y,Scale_Display=sc$Scale_Display[j],Scope=sc$Scope[j])
  cv=covars(sc$Is_Normed[j],TRUE)
  S13[[n]]=cbind(meta,Covariates=paste(cv,collapse=' + '),as.data.frame(as.list(corr(d,x,y,cv,paste('S13',x,y)))))
  cv=covars(sc$Is_Normed[j]);row=meta;row$Covariates=paste(cv,collapse=' + ')
  for(g in c('HS','non-HS')){
   z=corr(d[d$HS_status==g,],x,y,cv,paste('S15',g,x,y));z=as.list(z);names(z)=paste0(if(g=='HS')'HS_' else 'nonHS_',names(z));row=cbind(row,as.data.frame(z))
  }
  S15[[n]]=row
 }
 cat('Completed correlations:',i,'/',nrow(rois),'ROIs\n')
}
S13=bind_rows(S13)%>%group_by(Scale,Family)%>%mutate(Q_BH=p.adjust(P_perm,'BH'),FDR_tests=n())%>%ungroup()
S15=bind_rows(S15)%>%group_by(Scale,Family)%>%mutate(HS_Q_BH=p.adjust(HS_P_perm,'BH'),nonHS_Q_BH=p.adjust(nonHS_P_perm,'BH'),FDR_tests_per_group=n())%>%ungroup()
S14=list();n=0
for(s in segments)for(j in seq_len(nrow(sc))){
 n=n+1;x=paste0('ipsi_hp_',s);y=sc$Scale[j];cv=covars(sc$Is_Normed[j],TRUE)
 a=d[complete.cases(d[,c(x,y,cv,'ipsi_whole_hippocampus')]),];a$Y=as.numeric(a[[y]]);a$Segment=as.numeric(a[[x]])
 m0=lm(reformulate(c(cv,'ipsi_whole_hippocampus'),'Y'),a);m1=update(m0,.~.+Segment);an=anova(m0,m1);ss=summary(m1)
 # Diagnostic: residualized segment/whole collinearity after the other covariates.
 Q=qr(model.matrix(reformulate(cv),a));rc=cor(qr.resid(Q,a$Segment),qr.resid(Q,a$ipsi_whole_hippocampus))
 S14[[n]]=data.frame(Segment=x,Scale=y,Scale_Display=sc$Scale_Display[j],Scope=sc$Scope[j],n=nrow(a),Covariates=paste(cv,collapse=' + '),R2_reduced=summary(m0)$r.squared,R2_full=ss$r.squared,Delta_R2=ss$r.squared-summary(m0)$r.squared,F=an$F[2],df1=an$Df[2],df2=df.residual(m1),P_F=an$`Pr(>F)`[2],Delta_AIC=AIC(m1)-AIC(m0),Segment_beta=coef(m1)['Segment'],Segment_SE=ss$coefficients['Segment','Std. Error'],Segment_VIF=1/(1-rc^2))
}
# Explicit conservative correction over the complete 57-test nested-model set (19 core outcomes x 3 segments).
S14=bind_rows(S14)%>%mutate(Q_BH=p.adjust(P_F,'BH'),FDR_tests=n())
stopifnot(nrow(S13)==nrow(rois)*nrow(sc),nrow(S14)==length(segments)*nrow(sc),nrow(S15)==nrow(rois)*nrow(sc))  # 399 / 57 / 399 at 19 scales
write.csv(S13,file.path(rd,'Table_S13_Complete_HS_adjusted_correlations.csv'),row.names=FALSE,fileEncoding='UTF-8',na='')
write.csv(S14,file.path(rd,'Table_S14_Complete_segment_incremental_models.csv'),row.names=FALSE,fileEncoding='UTF-8',na='')
write.csv(S15,file.path(rd,'Table_S15_Complete_HS_nonHS_effect_sizes.csv'),row.names=FALSE,fileEncoding='UTF-8',na='')
write.csv(sc,file.path(vd,'S13_S15_Outcome_scope.csv'),row.names=FALSE,fileEncoding='UTF-8')
write.csv(rois,file.path(vd,'S13_S15_ROI_scope.csv'),row.names=FALSE,fileEncoding='UTF-8')
summary=bind_rows(S13%>%group_by(Family)%>%summarise(Table='S13',Tests=n(),Corrected=sum(Q_BH<.05),.groups='drop'),data.frame(Family='Ipsilateral hippocampal segments',Table='S14',Tests=nrow(S14),Corrected=sum(S14$Q_BH<.05)))
write.csv(summary,file.path(vd,'S13_S15_summary.csv'),row.names=FALSE,fileEncoding='UTF-8')
print(summary);print(S13%>%filter(ROI=='ipsi_whole_hippocampus',Q_BH<.05)%>%select(Scale_Display,n,r,P_perm,Q_BH))
print(S14%>%filter(Scale_Display=='WMS FSMQ')%>%select(Segment,n,Delta_R2,F,df2,P_F,Q_BH,Delta_AIC))
cat('S14 max/median delta R2:',max(S14$Delta_R2),median(S14$Delta_R2),' AIC increased:',sum(S14$Delta_AIC>0),'\n')
capture.output(sessionInfo(),file=file.path(vd,'S13_S15_R_sessionInfo.txt'))
