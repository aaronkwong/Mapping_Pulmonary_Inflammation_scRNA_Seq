#this function takes the path starting from within, the google drive folder and modifies the beggining of the path for the appropriate PC. For example
#to access the "Masters" folder, you would write: gdpath("Masters") 
#external defines whether the file specified is in the cloud folder. If external=TRUE, then the root is simply prepended to x and returned (without any other modification)
gdpath<-function(x,root,external=FALSE){
	if (x==""){
		stop("Error: no path specified.")
	}
	#prepend the root if provided
	if (missing(root)){
		prepend<-""
	}else{
		#if the root doesn't end in a slash, we should add one
		if(substr(root,nchar(root),nchar(root))!="/"){
			root<-paste0(root,"/")
		}
		prepend=root
	}

	if(external==TRUE){
		return(paste0(prepend,x,sep=""))
	}else{
		stop("This function cannot be run with external=FALSE")
	}
}

# #this function can souce libraries just by their name from the Rscript library
# #example of use: gdlib(c("win_slash","encryption_utilities"))
# #to list libraries available use gdlib("lib")
# gdlib<-function(x,script_library_path=gdpath("Masters/R Script Folder/lib"),script_library_database=gdpath("Masters/R Script Folder/Script_Manager/Script_Management_Data/Master_Script_Data.txt")){
# 	if(length(x)>1){
# 		for (package in x){
# 			source(paste0(script_library_path,"/",package,".R"))
# 		}
# 	}else if (length(x)==1 & x!="lib"){
# 		source(paste0(script_library_path,"/",x,".R"))
# 	}else if (x=="lib"){
# 		library_data<-read.delim(script_library_database,stringsAsFactors=F)
# 		print(library_data$Script_Name)
# 	}
# }
