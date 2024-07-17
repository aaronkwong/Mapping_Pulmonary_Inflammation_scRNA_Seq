#this function extracts the last object in a given path whether it is a directory or an object by grabbing everything after the last slash
last_obj_from_path<-function(x){
	start<-max(gregexpr("/",x)[[1]])+1
	end<-nchar(x)
	result<-substr(x,start,end)
	return(result)
}

path_without_last_directory<-function(x){
	start<-1
	end<-max(gregexpr("/",x)[[1]])
	result<-substr(x,start,end)
	return(result)
}

#this function removes the last object from a path (identical to "path_without_last_directory" function)
path_without_last_object<-function(x){
	start<-1
	end<-max(gregexpr("/",x)[[1]])
	result<-substr(x,start,end)
	return(result)
}
