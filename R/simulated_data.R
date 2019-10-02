#'Generate simulated data
#'@description It contains number of patients and risky population.
#'@docType data
#'@keywords datasets
#'@name simulated_data_1
#'@usage data(simulated_data_1)
#'@format  data frame
#'@examples
#'#Generate simulated data
#'#Given incidence rate function 1
#'FT1<-function(x,y){
#'  rate<-0.00148*sin(0.5*pi*x*y)+0.00002
#'  return(rate)
#'}
#'#Given risky population function
#'M<-function(x,y){
#'  r0=152040 ;r1=-285270; r2=110410; r3=173900
#'  r4=-49950 ;r5=-33630 ; r6=-19530; r7=-110330
#'  r8=88840  ;r9=-7990
#'  population<-r0+r1*x+r2*y+r3*x^2+
#'              r4*x*y+r5*y^2+r6*x^3+
#'              r7*x^2*y+r8*x*y^2+r9*y^3
#'  return(population)
#'}
#'ages<-35:85
#'years<-1988:2007
#'gen_data<-function(ages,years,FT,M){
#'   x<-mapping_to_01(ages)
#'   y<-mapping_to_01(years)
#'   disease<-outer(x,y,M)*outer(x,y,FT)
#'   population<-outer(x,y,M)
#'   zero<-matrix(0,ncol = 2*length(y),nrow = min(ages)-1)
#'   simulated_data<-rbind(zero,cbind(disease,population))
#'   colnames(simulated_data)<-rep(years,2)
#'   row.names(simulated_data)<-1:max(ages)
#'   return(simulated_data)
#'}
#'simulated_data_1<-gen_data(ages,years,FT1,M)
#'simulated_data_1
NULL
#'Generate simulated data
#'@description It contains number of patients and risky population.
#'@docType data
#'@keywords datasets
#'@name simulated_data_2
#'@usage data(simulated_data_2)
#'@format  data frame
#'@examples
#'#Generate simulated data
#'#Given incidence rate function 2
#'FT2<-function(x,y){
#'  rate<-0.00148*sin(0.5*pi*x*(y+0.2))+0.00002
#'  return(rate)
#'}
#'#Given population function
#'M<-function(x,y){
#'  r0=152040 ;r1=-285270; r2=110410; r3=173900
#'  r4=-49950 ;r5=-33630 ; r6=-19530; r7=-110330
#'  r8=88840  ;r9=-7990
#'  population<-r0+r1*x+r2*y+r3*x^2+
#'              r4*x*y+r5*y^2+r6*x^3+
#'              r7*x^2*y+r8*x*y^2+r9*y^3
#'  return(population)
#'}
#'ages<-35:85
#'years<-1988:2007
#'gen_data<-function(ages,years,FT,M){
#'   x<-mapping_to_01(ages)
#'   y<-mapping_to_01(years)
#'   disease<-outer(x,y,M)*outer(x,y,FT)
#'   population<-outer(x,y,M)
#'   zero<-matrix(0,ncol = 2*length(y),nrow = min(ages)-1)
#'   simulated_data<-rbind(zero,cbind(disease,population))
#'   colnames(simulated_data)<-rep(years,2)
#'   row.names(simulated_data)<-1:max(ages)
#'   return(simulated_data)
#'}
#'simulated_data_2<-gen_data(ages,years,FT2,M)
#'simulated_data_2
NULL
