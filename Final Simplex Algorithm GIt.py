#!/usr/bin/env python
# coding: utf-8

# # Matrix Entries

# # A

# In[ ]:


import numpy as np                                    #A-CONSTRAINT MATRIX
m=int(input("Eneter the number of rows"))
n=int(input("Enter the number of columns"))
A=[]
for i in range(m):
    a=[]
    for j in range(n):  #generating each row
        b=float(input("Enter the element(row={},col={}) of the constraint matrix".format(i+1,j+1)))
        a.append(b)
    A.append(a)         #attaching each row to the matrix A
A=np.asarray(A)
print(A)


# In[ ]:


x=int(input("Eneter the number of rows for H="))
H=[]
for i in range(x):
    h=[]
    for j in range(n):  #generating each row
        b=float(input("Enter the element(row={},col={}) of the constraint matrix".format(i+1,j+1)))
        h.append(b)
    H.append(h)         #attaching each row to the matrix H
H=np.asarray(H)
print(H)


# ### b

# In[ ]:


b=[]                                                   #B-SOLUTION MARTIX
for i in range(m):
    bc=[]
    z=float(input("Enter the elements of b="))
    bc.append(z)
    b.append(bc)
b=np.array(b)    
b=b.reshape((m,1))
print(b)


# In[ ]:


d=[]                                                   #B-SOLUTION MARTIX
for i in range(x):
    bc=[]
    z=float(input("Enter the elements of b"))
    bc.append(z)
    d.append(bc)
d=np.array(d)    
d=d.reshape((x,1))
print(d)


# # c

# In[ ]:


c=[]                                                 #OBJECTICVE FUNCTION
for i in range(n): #n elements in objective function
    h=float(input("Enter the coefficents of the objective function"))
    c.append(h)
c=np.asarray(c)        
print(c)


# #### bfs1

# In[ ]:


bsf=[]                                           #Input the first BFS
for i in range(n): #n elements in bsf
    s=float(input("Enter the element of BSF1"))
    bsf.append(s)


# In[ ]:


slack=[]              #adding slacks to A
for i in range(m):
    subslack=[]
    for j in range(m):
        if i==j:
            subslack.append(1)
        else:
            subslack.append(0)
    slack.append(subslack) 
slack=np.asarray(slack)
print(A,slack)
A=np.append(A, slack, axis=1)
H=np.array(H)       #H integration with A
dh=H.shape
z=dh[0]
duh=[]
for i in range(z):
    forh=[]
    for j in range(m): #old 'm' should be used here
        forh.append(0)
    duh.append(forh) 
duh=np.asarray(duh)    
H=np.append(H,duh,axis=1)
H=np.asarray(H)
A=np.asarray(A)
A=np.append(A,H,axis=0)
print(A)
b=np.append(b,d,axis=0)
d=A.shape
m1=d[0]
n1=d[1]
    


# In[ ]:


def gaussianeleminaiton(A):
    z=max(m,n)
    for i in range(z):
        x=A[i,i]
        if x==0:                         # shuffling rows to get pivot on top.
            for j in range(i+1,m):
                if A[j,i]!=0:
                    A=list(A)
                    A[i],A[j]=A[j],A[i]
                    A=np.asarray(A)
                    x=A[i,i]
        for col in range(n):             #First operation 1st element to 1
            A[i,col]=A[i,col]/x
        for row in range(i+1,m):           #[1,0]=0
            y=A[row,i]/A[i,i]
            for col in range(n):
                A[row,col]=A[row,col]-(y*A[i,col])
            
d=A.shape
m1=d[0]
n1=d[1]


# In[ ]:


print(m1,n1)
slacka=[]              #adding slacks to A
for i in range(m1):
    subslack=[]
    for j in range(m1):
        if i==j:
            subslack.append(1)
        else:
            subslack.append(0)
    slacka.append(subslack) 
slacka=np.asarray(slacka)
A=np.append(A, slacka, axis=1)
d=A.shape
m=d[0]
n=d[1]
print(m)


# In[ ]:


print(A,b)


# In[ ]:


for i in range(m):   #ensuring positivity of auxiallary lp
    if b[i]<0:
        b[i]=((-1)*b[i])
        for j in range(n):
            A[i,j]=((-1)*A[i,j])
print(A,b,m,n)            


# In[ ]:


print(m,m1)
absf=[]
for i in range(n):
    if i<=n1-1:
        absf.append(0)
    else:
        absf.append(1)
print(abfs)
absf=np.asarray(absf)


# In[ ]:


c=[]
for i in range(n):
    if i<=n1-1:
        c.append(0)
    else:
        c.append(1)
print(c)
c=np.asarray(c)


# # Pivot Max coefficient

# In[ ]:


def pivotmax(A,b,c,bsf):
    cb=[]                                          #Reduced cost coefficient elements
    cn=[]
    B=[]
    N=[]
    crindex=[]
    crindexb=[]
    count=0
    counter=0
    for i in range(n): #getting the indices of the basic and non basic variables
        if bsf[i]==0:
            crindex.append(i)  #indexes of non basic variables from A matrix
            cn.append(c[i])   
            N.append(A[:,i])
            count=count+1
        else:
            crindexb.append(i)
            cb.append(c[i])
            B.append(A[:,i])
            counter=counter+1
    cb=np.asarray(cb) #cb
    cn=np.asarray(cn) #cn
    B=np.asarray(B)  #B
    B=B.transpose()
    N=np.asarray(N) #N
    N=N.transpose()
    
    Binv=np.linalg.inv(B)#Binverse       #developing the reduced cost square coefficeint
    cbBi=np.matmul(cb, Binv)
    cr2=np.matmul(cbBi, N)#cb*Binverse*N
    cr=np.subtract(cn,cr2) #FINAL CR
    
    counter=0
    for i in cr:                        #OPTIMALITY CHECK
        if i<=0:
            counter=counter+1
    if counter==n-m:       #number of elements in cr = n-m
        ofv1=np.matmul(cb,Binv)
        ofv=np.matmul(ofv1,b)
        print("this is the optimal bsf, the bsf is ", bsf, "the objective funstion value is", ofv)
        x="optimal"
        return x
    
        
    check=[]                    #UNBOUNDEDNESS CHECK
    for i in range(n-m):       #range=n-m beacuse of the number of elements in RCCV 
        if cr[i]>=0:
            kr=N[:,i].reshape(m,1)
            count=0
            unbv=np.matmul(Binv,kr)#Binv*Kr is a m*1 matrix
            for i in unbv:
                if i<=0:
                    count=count+1
                    if count==m: #checking if all elements in unbv are negetive
                        check.append(1)
                
    if len(check)>0:    #checking if atleast one of Binv*Kr is negetive                           
        x="unbounded"
        return x  
        
        
    crl=list(cr)
    q=max(crl)      
    o=crl.index(q)    #index of max cr
    v=crindex[o]#finding the index (v) of the variable that will be turning from non basic to basic 
    shiftn=A[:,v]
    shiftn=np.reshape(shiftn,(m,1)) #the vector being used while shifting
      
    
    Binvb=np.matmul(Binv,b)
    Binvshiftn=np.matmul(Binv, shiftn)
    newxb=np.subtract(Binvb,Binvshiftn)
    
    
    comp=list(Binvb-newxb) #selecting the new bsf
    newxbl=list(newxb)
    Binvbl=list(Binvb)
    for i in range(m):
        if comp[i]>0 and newxbl[i]==min(newxbl):
            h=i  #index of the basic variable that will be turnig non basic.
            break
        else:
            continue
    
    j=crindexb[h]
    bsf=list(bsf)
    bsf[v]=1
    bsf[j]=0
    print(cr)
    print(bsf)
    pivotmax(A,b,c,bsf)
    
            
            

    
    


# # Pivot Bland

# In[ ]:


def pivotbland(A,b,c,bsf):
    cb=[]                                          #Reduced cost coefficient elements
    cn=[]
    B=[]
    N=[]
    crindex=[]
    crindexb=[]
    count=0
    counter=0
    for i in range(n): #getting the indices of the basic and non basic variables
        if bsf[i]==0:
            crindex.append(i)  #indexes of non basic variables from A matrix
            cn.append(c[i])   
            N.append(A[:,i])
            count=count+1
        else:
            crindexb.append(i)
            cb.append(c[i])
            B.append(A[:,i])
            counter=counter+1
    cb=np.asarray(cb) #cb
    cn=np.asarray(cn) #cn
    B=np.asarray(B)  #B
    B=B.transpose()
    N=np.asarray(N) #N
    N=N.transpose()
    
    Binv=np.linalg.inv(B)#Binverse       #developing the reduced cost square coefficeint
    cbBi=np.matmul(cb, Binv)
    cr2=np.matmul(cbBi, N)#cb*Binverse*N
    cr=np.subtract(cn,cr2) #FINAL CR
    
    counter=0
    for i in cr:                        #OPTIMALITY CHECK
        if i<=0:
            counter=counter+1
    if counter==n-m:       #number of elements in cr = n-m
        ofv1=np.matmul(cb,Binv)
        ofv=np.matmul(ofv1,b)
        print("this is the optimal bsf, the bsf is ", bsf, "the objective funstion value is", ofv)
        x="optimal"
        return x
    
        
    check=[]                    #UNBOUNDEDNESS CHECK
    for i in range(n-m):       #range=n-m beacuse of the number of elements in RCCV 
        if cr[i]>=0:
            kr=N[:,i].reshape(m,1)
            print("kr",kr)
            count=0
            unbv=np.matmul(Binv,kr)#Binv*Kr is a m*1 matrix
            for i in unbv:
                if i<=0:
                    count=count+1
                    if count==m: #checking if all elements in unbv are negetive
                        check.append(1)
                
    if len(check)>0:    #checking if atleast one of Binv*Kr is negetive                           
        x="unbounded"
        return x  
        
        
    crl=list(cr)
    for i in range(n-m):
        if crl[i]>0:
            o=i
            break
    v=crindex[o]#finding the index (v) of the variable that will be turning from non basic to basic 
    shiftn=A[:,v]
    shiftn=np.reshape(shiftn,(m,1)) #the vector being used while shifting
      
    
    Binvb=np.matmul(Binv,b)
    Binvshiftn=np.matmul(Binv, shiftn)
    newxb=np.subtract(Binvb,Binvshiftn)
    print("Binvb", Binvb)
    print("newxb", newxb)
    
    
    comp=list(Binvb-newxb) #selecting the new bsf
    newxbl=list(newxb)
    Binvbl=list(Binvb)
    for i in range(m):
        if comp[i]>0 and newxbl[i]==min(newxbl):
            print(comp[i])
            h=i  #index of the basic variable that will be turnig non basic.
            break
        else:
            continue
    
    j=crindexb[h]
    bsf=list(bsf)
    bsf[v]=1
    bsf[j]=0
    print(cr)
    print(bsf)
    pivotbland(A,b,c,bsf)
    
            
            

    
    


# # functions run here

# In[ ]:


pivotbland(A,b,c,bsf)


# In[ ]:


pivotmax(A,b,c,bsf)


# In[ ]:





# In[ ]:




