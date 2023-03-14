function Z = Maximise12(A,b,C,Opt);
Q=size(A,1);
P=size(A,2);%Loads the size of A in both dimensions into variables


M=eye(Q);
k=1;
V=[];
IndependentVars=[];
for i=1:P %This checks our A to see which, if any, artificial variables to add. We do this by finding any variables that appear in only a single equation and excluding columns of the QxQ matrix that have a 1 in the same row.
    if sum(A(:,i))==1
        OnesInColumn=0;
        badcolumn=0;
        for j=1:Q;
           if A(j,i)==1
               OnesInColumn=OnesInColumn+1;
               RowIndex=j;
           elseif A(j,i)~=0
               badcolumn=1;
               break
           end
        end
        if OnesInColumn==1 && badcolumn==0
                V(k)=RowIndex; %Storing which row the 1 is in 
                IndependentVars(k)=i; %Storing which column the 1 is in
                k=k+1;
        end
    end
end
M(:,V)=[];
L=size(M,2);
ArtificialVars=[P+1:P+L]; %This is stored so that one can check at the end if the artificial variables are nonbasic, if they are basic then the problem should have no feasible solutions.
A(:,P+1:P+L)=M
if Opt=="Min"
    C=-C;
    O=-1;
elseif Opt=="Max"
    O=1;
end
%Our code only solves the maximisation problem so this bit above converts a
%minimisation problem into a maximisation problem in order to proceed. Any
%further use of the variable O is being used for this purpose also.

C(P+1:P+L,:)=-1000000*ones(L,1);
P=size(A,2);
W=P-Q+1;
% P is the number of columns of A and Q is the number of rows of A
X=IndependentVars; %Our vector X keeps track of which variables are basic and their order in the tableau
if length(X)~=Q %This appends any necessary artifical variables to X
    X((length(X)+1):Q)=(P-(Q-length(X))+1):P;
end
B=A(:,X);
% This is our starting B by constructing a basic starting solution,
% specifically the basic starting solution requiring the minimum number of
% artificial variables.
CB=C(X);
% This is our new C value based on the variables in B
XB=B\b;
Zrow=(CB'*(B\A)-C');
% Zrow is our Z row based on the initial variables
Z=dot(CB,XB);
% Calculates our Z value based on the chosen variables
Zrow(Zrow>=0)=nan; %Replacing the positive elements of the Z row so that we can find the most negative number
MinZrow=min(Zrow) %Finds the largest magnitude negative number in our Z row
X
while MinZrow<0 % Checking our condition for optimisation
    [~,tempenterindex]=min(Zrow) %Stores which variable xi is the most negative in the Z row
    enterindexvector=1:P;
    enterindexvector(X)=[]
    enterindex=enterindexvector(tempenterindex)
    SolutionCheck=B\A(:,enterindex)
    
    r1=0; %r1 and r2 are placeholder variables that 
    Check=[];
    Checkindex=[];
    for i = 1:Q
        if SolutionCheck(i)>0
            r1=r1+1;
            Check(r1)=XB(i)/SolutionCheck(i)
            Checkindex(r1)=i
        end
    end
    
    r2=0;
    PosCheck=[];
    PosCheckindex=[];
    for i = 1:r1
        if Check(i)>=0
            r2=r2+1;
            PosCheck(r2)=Check(i)
            PosCheckindex(r2)=Checkindex(i)
        end
    end
    [~,templeaveindex]=min(PosCheck);
    leaveindex=PosCheckindex(templeaveindex)
    if length(leaveindex)==0
        break
    end
    %We now reconstruct all of our variables by swapping out the leaving variable for the entering variable
    B(:,leaveindex)=A(:,enterindex)
    X(leaveindex)=enterindex
    CB=C(X);
    XB=B\b;
    Zrow=(CB'*(B\A)-C');
    Z=O*dot(CB,XB)
    Zrow(Zrow>=0)=NaN;
    Zrow=Zrow;
    Zrow(X)=[]
    MinZrow=min(Zrow) % Now we have completed the first iteration and if it isn't optimal, the while loop will repeat            
end

z2=XB.*CB;
ArtificialCheck=intersect(ArtificialVars,X);
if (isnan(Z) || isinf(Z)) && isempty(ArtificialCheck)
    "The problem has no finite solutions."
end
if ~isempty(ArtificialCheck) %This tells us if an artificial variable is still present in our solution. If there is then we know that there was no feasibe solution to the problem
    "The problem has no feasible solutions."
end
if (~isnan(Z) && ~isinf(Z)) && isempty(ArtificialCheck) %If there is a finite feasible solution, this causes our function to output: which variables are still being used, the values they take, and the optimum Z value that comes from them
    XB
    X
    Z
end
end