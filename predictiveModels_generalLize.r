
### 
rm(list=ls());
library(gbm) #generalized boosted models
library(randomForest) #random forest
library(ROCR) #for ROC curves
library(ggplot2) #for nice graphs
library(openxlsx);
library(dplyr);
library(reshape2);
library(caret);
library(impute);
library(ggpubr);



MyDir <- "H:/tmp/ML"
setwd(MyDir)


### processing inhouse data
## gene-protein mapping through DAVID
dd_exp=read.xlsx("Source Data.xlsx", sheet="expression_164proteins", startRow=1, check.names=F, sep.names="_");
genes_164=dd_exp[, c("Protein_Accession", "geneID", "Gene_Symbol")];

### check the values in the first column are unique or not to decide if need aggregate at gene levels
length(as.character(dd_exp[, 1]))==length(unique(as.character(dd_exp[, 1])))
### convert protein level at gene levels
dd_exp = dd_exp[, c(-1, -2)];
### data transport 
dd_exp_data=dd_exp[, -1];
rownames(dd_exp_data)=as.character(dd_exp[, 1]);
exp_data=t(dd_exp_data);
colnames(exp_data)

### processing CPTAC data
### check how many protein-coding genes of 164 proteins were detected in CPTAC cohort
## generate gene protein mapping through DAVID in CPTAC
gene_mapping_CPTAC =read.xlsx("Source Data.xlsx", sheet="protein-gene_mapping_CPTAC",startRow=1, check.names=F, sep.names="_");

gene_mapping_164_CPTAC=merge(gene_mapping_CPTAC, genes_164, by.x="To", by.y="geneID");
print(paste("number of genes in CPTAC: ",  length(unique(gene_mapping_164_CPTAC[, "Gene_Symbol"])), sep=""))
gene_mapping_164_CPTAC=gene_mapping_164_CPTAC[, c("From", "Gene_Symbol")];

unique_genes=unique(gene_mapping_164_CPTAC[, "Gene_Symbol"]);

dd_CPTAC=read.xlsx("Source Data.xlsx", startRow=1, sheet="Global-Proteome-G3_CPTAC", check.names=F, sep.names="_");
head(dd_CPTAC)
dd_CPTAC = dd_CPTAC[, c(1, 13:ncol(dd_CPTAC))];
dd_CPTAC=merge(gene_mapping_164_CPTAC, dd_CPTAC, by.x=1, by.y=1);

### aggregate protein expression data at gene level
dd_CPTAC_agg=aggregate(dd_CPTAC[, -c(1, 2)], by=list(dd_CPTAC[, "Gene_Symbol"]), median, na.rm=TRUE);
colnames(dd_CPTAC_agg)[1]="Genes";

### convert sample names 
which_TCGA=grep("TCGA", colnames(dd_CPTAC_agg));
dd_CPTAC_agg=dd_CPTAC_agg[, c(1, which_TCGA)];
colnamess=colnames(dd_CPTAC_agg);
name_info=colsplit(colnamess, pattern="\\.", names=c("T1", "T2"));
namess=paste("TCGA-", name_info[, 1], sep="");
namess=c("Genes", namess[-1]);
colnames(dd_CPTAC_agg)=namess;

### aggragate duplicates expression at case level
which_dup=which(duplicated(namess));
sample_dup=unique(namess[which_dup]);

#### find samples without duplicates
which_Dups=which(namess %in% sample_dup)
dd_CPTAC_noDup=dd_CPTAC_agg[, -which_Dups]; ## include Genes column

### merge data for duplicates
samples_dup=c();
for (i in 1:length(sample_dup)){
   sample_name = sample_dup[i];
   dd_dup = dd_CPTAC_agg[, which(namess == sample_name)];
   dd_one=apply(dd_dup, 1, median, na.rm=TRUE);
   samples_dup=cbind(samples_dup, dd_one);
}

colnames(samples_dup)=sample_dup;

### generate normalized data at gene level per case
dd_CPTAC_data = cbind(dd_CPTAC_noDup, samples_dup);

print(dim(dd_CPTAC_data))

dd_CPTAC_data[, -1] = apply(dd_CPTAC_data[, -1], 2, as.numeric);
rownames(dd_CPTAC_data)=as.character(dd_CPTAC_data[, 1]);
CPTAC_data=t(dd_CPTAC_data[, -1]);
genes=colnames(CPTAC_data);


samples_CPTAC=rownames(CPTAC_data);
CPTAC_data=cbind(samples_CPTAC, CPTAC_data);

### get the IHC-based subtypes
dd_clinical_BRCA = read.xlsx("Source Data.xlsx", sheet="TCGA_clinical_subtypes", startRow=1, check.names=F, sep.names="_");

### extract samples with Luminal and TN subtype
dd_clinical_BRCA= dd_clinical_BRCA[, c("samples", "Luminal-TN")] %>% dplyr::filter(`Luminal-TN` %in% c("Luminal", "TN"));
CPTAC_data_subtype =merge(dd_clinical_BRCA, CPTAC_data, by.x=1, by.y=1);

### convert into 0 and 1 for factor
CPTAC_data_subtype[, "subtype"]=1*(CPTAC_data_subtype[, "Luminal-TN"] == "TN");

### re-arange the data with subtype as first column and genes as the follow columns
CPTAC_data_subtypes=data.frame(CPTAC_data_subtype[, c("subtype", genes)], check.names=F);
rownames(CPTAC_data_subtypes)=as.character(CPTAC_data_subtype[, "samples"]);
CPTAC_data_subtypes[, -1]=apply(CPTAC_data_subtypes[, -1], 2, as.numeric);

### check the missing value frequency per gene: which genes expressed at 80% samples
number_exp=apply(CPTAC_data_subtypes[, -1], 2, function(x) {return(length(which(!is.na(x))));})
which_80perc=which(number_exp/nrow(CPTAC_data_subtypes)>0.8);

length(which_80perc)==(ncol(CPTAC_data_subtypes)-1) ### if equal, no need to remove genes

samples=rownames(CPTAC_data_subtypes);
subtype=CPTAC_data_subtypes[, 1];
d1=cbind(samples, subtype);

d2=t(CPTAC_data_subtypes[, -1]);

### impute genes for missing value
d2_imputed=impute.knn(d2, k=10, rowmax=0.8, colmax=0.8, maxp=nrow(d2), rng.seed=1)$data;

d22=t(d2_imputed)
samples=rownames(d22);
d22=cbind(samples, d22);
d12=merge(d1, d22, by.x=1, by.y=1);


CPTAC_final=data.frame(apply(d12[, -1], 2, as.numeric));
CPTAC_final=CPTAC_final[, c("subtype", genes)];


### re-arranage the in-house data
exp_data = exp_data[, c("cohort", "Luminal-TN", genes)];
exp_data=data.frame(exp_data, check.names=F);
exp_data[, -c(1, 2)]=apply(exp_data[, -c(1, 2)], 2, as.numeric);

exp_data = exp_data %>% dplyr::filter(cohort %in% c("Training", "Testing"));

### check the levels of subtype
table(exp_data[, "Luminal-TN"])
exp_data[, "Luminal-TN"] = factor(exp_data[, "Luminal-TN"], levels=c("Luminal", "TN"), ordered=TRUE);
#for factor variables, encode them as 0/1 continuous variables

exp_data[, "subtype"]=1*(exp_data[, "Luminal-TN"] == "TN");

### check
table(exp_data[, c("subtype", "Luminal-TN")])

### generate training and testing data 
dd_exp_training=exp_data[which(exp_data[, "cohort"] == "Training"), c("subtype", genes)];
dd_exp_testing=exp_data[which(exp_data[, "cohort"] == "Testing"), c("subtype", genes)];

#Begin to build models
set.seed(1)
#Logistic regression
logitMod1 <- glm(formula=subtype ~ ., data=dd_exp_training, family="binomial")
summary(logitMod1)

#Random Forest
rfMod2 <- randomForest(formula=as.factor(subtype) ~ ., data=dd_exp_training, ntree=1000, importance=TRUE) #need as factor, else uses linear regression
print(rfMod2)


varImpPlot(rfMod2, type=1) #importance plot, the Scores with the strongest predictors
					   
#generalized boosted regression
gbMod3 <- gbm(formula=subtype ~ ., distribution="bernoulli", data=dd_exp_training, interaction.depth=5)
summary(gbMod3) # ANP32E is the main factor, other variables contribute much less

#Apply model to the testing data set to determine the model and Predictions on the testing data set
## Apply the determined model to validation cohorts, generate ROC curves

#predicted probabilities
new_datasets=c("Training", "Testing", "CPTAC");

plots=list();
for (i in 1:length(new_datasets)){
    new_db=new_datasets[i];
    
	if (new_db == "Training") {
	   new_dataset=dd_exp_training;
	} 
	
	if (new_db == "Testing") {
	   new_dataset=dd_exp_testing;
	} 
	
	if (new_db == "CPTAC") {
	  new_dataset = CPTAC_final;
	}
	PredLogit <- predict(logitMod1, newdata=new_dataset, type="response")
	PredRF <- predict(rfMod2, newdata=new_dataset, type="prob")[,2] #response for this produces 0/1
	PredGB <- predict(gbMod3, newdata=new_dataset, type="response", n.trees=100)

	#  predict at >5% defined as TN, otherwise as Luminal
	#rows are observed, columns are predicted
	t1=table(new_dataset$subtype,(PredLogit > 0.05)*1) #Logit, does not work on this, because overfitting in training dataset
	t2=table(new_dataset$subtype,(PredRF    > 0.05)*1) #Random Forest
	t3=table(new_dataset$subtype,(PredGB    > 0.05)*1)

	#hist(PredRF) 
	#hist(PredGB)

	#### generate ROC curves for each of these predictors
	PredLogit_P <- prediction(PredLogit, new_dataset$subtype)
	PredRF_P <- prediction(PredRF, new_dataset$subtype)
	PredGB_P <- prediction(PredGB, new_dataset$subtype)

	#performance of the ROC curves
	L1 <- performance(PredLogit_P, "tpr", "fpr")
	R2 <- performance(PredRF_P, "tpr", "fpr")
	G3 <- performance(PredGB_P, "tpr", "fpr")

	#convert into a dataframe
	ROC_data <- data.frame(rbind(cbind(L1@x.values[[1]],L1@y.values[[1]],1),
                             cbind(R2@x.values[[1]],R2@y.values[[1]],2),
                             cbind(G3@x.values[[1]],G3@y.values[[1]],3))
	        )
	names(ROC_data) <- c("fpr","tpr","Model")
	ROC_data$Model <- as.factor(ROC_data$Model)
	levels(ROC_data$Model) <- c("Logit","Random Forest","GBM")

	#area under the curve for each method
	auc_Logit = performance(PredLogit_P, "auc")@y.values[[1]]
	auc_RF = performance(PredRF_P, "auc")@y.values[[1]]
	auc_GB = performance(PredGB_P, "auc")@y.values[[1]]

	note1=paste("AUC(Logit)=", signif(auc_Logit, digits=2), sep="");
	note2=paste("AUC(RF)=", signif(auc_RF, digits=2), sep="");
	note3=paste("AUC(GB)=", signif(auc_GB, digits=2), sep="");


	colorss = c("Logit" = "royalblue", "Random Forest" = "orange", "GBM" = "hotpink");
	p <- ggplot(data=ROC_data, aes(x=fpr,y=tpr,color=Model)) + geom_abline(slope=1, size=0.5, color="gray") + geom_line(size=0.5) 
	p = p + scale_color_manual(values = colorss)
	p = p +xlab("False Positive Rate") + ylab("True Positive Rate") + ggtitle(new_db)
	p=p+ theme(plot.title = element_text(size=7),
	          axis.title = element_text(size = 6), 
			  axis.text.x = element_text(size = 6), 
			  axis.text.y = element_text(size = 6), 
			  legend.text = element_text(size = 6), 
			  legend.title = element_text(size = 7))
    # p= p+ theme(text = element_text(size = 1))
	p = p + geom_text(x=0.85, y=0.25, label=note1, size=2,  col="royalblue")
	p = p + geom_text(x=0.85, y=0.15, label=note2, size=2,  col="orange")
	p = p + geom_text(x=0.85, y=0.05, label=note3, size=2, col="hotpink")
	plots[[i]]=p
    png(paste("ROC_", new_db, ".png", sep=""), height=2, width=4, units="in", res=300);
	print(p);
	dev.off();
}


p1=plots[[1]];
p2=plots[[2]];
p3=plots[[3]];
ggg=ggarrange(p1, p2, p3,  ncol = 3,  nrow=1, labels=c("a", "b", "c"),  common.legend =TRUE)
png("ROC_training_testing_validation.png", height=3, width=7, units="in", res=300);
print(ggg);
dev.off();
dev.off();
	
