# update the metadata field of a given object
SC.UpdateMeta<-function(meta, categories, keys, values) {

keys<-keys[!is.na(keys)];
keys<-keys[keys!=''];

if (length(keys)<1) cat('key is missing, no metadata updates\n')
else if (length(categories)<1) cat('categories is missing, no metadata updating\n')
else if (length(values)<1) cate('value is missincategoriesg, no metadata updating\n')
else {
if(length(values)<length(keys)) value[(length(values)+1):length(keys)]<-NA;
if(length(categories)<length(keys)) categories[(length(categories)+1):length(keys)]<-categories[length(categories)];

for (i in 1:length(keys)) {
if (is.null(meta[[categories[i]]])) { # the category not exists, create new
x<-values[i];
names(x)<-keys[i];
meta[[length(meta)+1]]<-x;
names(meta)[length(meta)]<-categories[i];
}
else { # update key-value pair in existing categories
meta[[categories[i]]][keys[i]]<-values[i];
}
} # end of for
meta;

}

}

