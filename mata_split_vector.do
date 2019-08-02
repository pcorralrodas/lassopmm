clear all
set more off

mata


y=rnormal(103,1,0,1)

function _br(y){
	pointer(real matrix) rowvector kk
	kk = J(1,10,NULL)
 
 add=floor(rows(y)/10)

 final=mod(rows(y),10)
 
 kk=J(1,10,NULL)
 a=1
 
	for(i=1;i<=rows(y);i++){
		nexti = i+add-1
		kk[a] = &(y[|i,1\nexti,1|])

		if (a==10){
			kk[a] = &(y[i..(nexti+final)])			
		}
		a=a+1
		i=nexti
	}
 
}

G=_br(y)