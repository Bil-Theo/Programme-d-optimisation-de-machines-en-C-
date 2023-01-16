#include <iostream>
#include <list>
#include <list>
#include <iomanip>
using namespace std; 


int it; //nombre d'etirations
static int max(int a, int b){
	int c=a;
	if(b>a)c=b;
	return c;
}

int getElement(list <int> l, int position){
    int val;
    if(l.size()<position){//si la position demander ne contient aucun élement 
        cout<<"l'element à la position "<<position<<"n'existe pas"<<endl;
        //arrêter le programme pour éviter les problèmes d'execution
        exit(1);
    }
    else{
        list <int> temp;
        //on copie la liste l dans temp pour ne pas perdre l
        temp = l;
        //on supprime les élement de la list jusqu'à l'evant qui precède position
        for(int i = 0; i<position; i++) temp.pop_front();
        //on recupère la valeur de la position dans val on retourne ca
        val = temp.front();
    }
    return val;
    
}

int getElement(list<list<int>> l, int posLine, int posCol){
    int val;
    if(l.size()<posLine){
        //on verifie si la position demander de ligne existe
         cout<<"la ligne "<<posLine<<"n'existe pas"<<endl;
        //arrêter le programme pour éviter les problèmes d'execution
        exit(1);
    }
    else if(l.front().size()<posCol){
        //on verifie si la position demander de colonne existe
         cout<<"la colone "<<posCol<<"n'existe pas"<<endl;
        //arrêter le programme pour éviter les problèmes d'execution
        exit(1);
    }
    else{
        list<list<int>> t;
        //on recopie la matrice de list l dans t pour ne pas perdre l
        t = l ;
        //on supprime les ligne de t (copie de l) jusqu'à arriver à la ligne demander pour extraire
        for(int i = 0; i<posLine; i++) t.pop_front();
        //on supprime les elements(colonnes) de la lignes demander jusqu'à arriver à l'element demander
        for(int i = 0; i<posCol; i++) t.front().pop_front();
        //on extrait l'element demander on met dans val puis on le retourne
        val = t.front().front();
    }
    
    return val;
    
}

void setElement(list<int> &l, int position, int elt){
    
    list<int> temp;
    temp = l;
    int var_help;
    if(l.size()<position){
        //on verifie si l'element qu'on veut modifier existe
        cout<<"L'element demander à la position "<<position<<"n'existe pas"<<endl;
        //on retourne l'ancienne liste sans modification
        return;
    }
    else{
        //on supprime les elements qui se trouve avant et la position demandée
        for(int i = 0; i<=position; i++) temp.pop_front();
        
        //on ajoute l'element à modifier à la place de la position demandée
        temp.push_front(elt);
        
        //on remet les elements supprimée avant pour avoir de nouveau la liste
        for(int i = 0; i<position; i++){
            var_help = getElement(l,i);
            temp.push_front(var_help);
        }
        
        
        //on retourn la liste modifié
        l = temp;
        return;
    }
}

void setElement(list<list<int>> &l, int posLine, int posCol, int elt){
    
    if(l.size()<posLine){
        cout<<"la ligne demandée "<<posLine<<" pour modifier l'element n'existe pas"<<endl;
        return;
    }
    else if(l.front().size()<posCol){
        
        cout<<"l'element demandée pour modifier l'element n'existe pas"<<endl;
        return;
    }
    else{
        list<list<int>> t;
        list<int> h;
        t = l;
        int var_help;
        
        for(int i = 0; i<posLine; i++) t.pop_front();
        h = t.front();
        for(int j = 0; j<=posCol; j++) t.front().pop_front();
        
        t.front().push_front(elt);
        
        for(int j = 0; j<posCol; j++){
            var_help = getElement(h,j);
            t.front().push_front(var_help);
        }
        
        for(int i = 0; i<posLine; i++){
            t.push_front(l.front());
            l.pop_front();
        }
        
        l = t;
        return;
    }
}

void saisir (int n, list <int> &vec){
	int e;
	for(int i=0; i<n; i++){
		cin>> e;
		vec.push_back(e);
	}
}

void saisir(int n, int m, list <list<int> > &p){
   list<int> vect;
	int e;
	for(int j=0; j<m; j++){
		cout<<"les temps de traitement sur M"<<j+1<<endl;
		vect.clear();
		for(int i=0; i<n; i++){
			cin>>e;
			vect.push_back(e);
		}
		p.push_back(vect);
	}
}

void affichage (int n, int tab[]){
	for(int i=0; i<n; i++)
		cout << tab[i] << " ";
		cout<< " "<<endl;
}

void permutation (int &a, int &b){
	int val=a;
	    a=b;
	    b=val;
}

void tri(int n, int tab[], int ind[]){
	int t[n];
	for(int i=0;i<n;i++) {
		ind[i]=i;
		t[i]=tab[i];
	}
	for(int i=0; i<n; i++){
		for(int j=i+1; j<n; j++){
			if(t[i]<t[j]){
				permutation(t[i],t[j]);
				permutation(ind[i],ind[j]);
			}	
		}
	}
}

void Tri(int n, int tab[], int ind[]){
	int t[n];
	for(int i=0; i<n; i++)	{
		//ind[i]=i;
		t[i]=tab[i];
	}
	for(int i=0; i<n; i++){
		for(int j=i+1; j<n; j++){
			if(t[i]<t[j]){
				permutation(t[i],t[j]);
				permutation(ind[i],ind[j]);
			}	
		}
	}
}

int ordonnancer (int n,  int n1, int m, list <list<int> > p, list<int> ind){
	int Cmax;
	int c[m][n1];
	int a, b, x, d, y;
	//initialiser la matrice "c" a 0
	for (int j=0;j<m;j++){ for(int i=0;i<n1;i++)	c[j][i]=0;	}
	//premiere tache
	for (int j=0;j<m;j++){
		for(int k=0;k<=j;k++){
			//c[j][ind[0]] = c[j][ind[0]] + p[k][ind[0]];
			a = getElement(ind, 0);
			b =  getElement(p, k, a);
			c[j][a] = c[j][a] + b;
		}
	}
	//premiere machine 
	//c[0][ind[0]]=0;
	c[0][a] = 0;
	for (int i=0;i<n;i++){
		for(int k=0;k<=i;k++){
			//c[0][ind[i]] = c[0][ind[i]] + p[0][ind[k]];
			a = getElement(ind, i);
			x = getElement(ind, k);
			b =  getElement(p, 0, x);
			c[0][a] = c[0][a] + b;
		}
	}
	for(int i=1; i<n;i++){
		for(int j=1; j<m;j++){
			//c[j][ind[i]] = max(c[j][ind[i-1]], c[j-1][ind[i]]) + p[j][ind[i]];
			a = getElement(ind, i);
			b =  getElement(p, j, a);
			d = getElement(ind, i-1);
			c[j][a] = max(c[j][d], c[j-1][a]) + b;

		}
	}
	//Cmax=c[m-1][ind[n-1]];
	x = getElement(ind, n-1);
	Cmax = c[m-1][x];
	//cout<<"Cmax= "<<Cmax<<endl;
	return Cmax;
}

int ordonnancer (int n,  int m, list <list<int> > p, list<int> ind){
	int Cmax;
	int c[m][n];
	int a,b,x,d,y;
	//initialiser la matrice "c" a 0
	for (int j=0;j<m;j++){ for(int i=0;i<n;i++)	c[j][i]=0;	}
	//premiere tache
	for (int j=0;j<m;j++){
		for(int k=0;k<=j;k++){
			//c[j][ind[0]] = c[j][ind[0]] + p[k][ind[0]];
			a = getElement(ind, 0);
			b =  getElement(p, k, a);
			c[j][a] = c[j][a] + b;
		}
		x = getElement(ind, 0);
		cout<<"T"<<x+1<<" sur M"<<j+1<<" : ["<<c[j][x]-getElement(p, j, a)<<" , "<<c[j][a]<<" ]"<<endl;
	}
	//premiere machine 
	x = getElement(ind, 0);
	c[0][x]=0;
	c[0][a] = 0;
	for (int i=0;i<n;i++){
		for(int k=0;k<=i;k++){
			//c[0][ind[i]] = c[0][ind[i]] + p[0][ind[k]];
			a = getElement(ind, i);
			x = getElement(ind, k);
			b =  getElement(p, 0, x);
			c[0][a] = c[0][a] + b;
		}
		y = getElement(ind, i);
		cout<<"T"<<y+1<<" sur M1  : ["<<c[0][y]-getElement(p, 0, y)<<" , "<<c[0][y]<<" ]"<<endl;
	}
	for(int i=1; i<n;i++){
		for(int j=1; j<m;j++){
			//c[j][ind[i]] = max(c[j][ind[i-1]], c[j-1][ind[i]]) + p[j][ind[i]];
			a = getElement(ind, i);
			b =  getElement(p, j, a);
			d = getElement(ind, i-1);
			c[j][a] = max(c[j][d], c[j-1][a]) + b;
			y = getElement(ind, i);
			cout<<"T"<<y+1<<" sur M"<<j+1<<" : ["<<c[j][y]-getElement(p, j, y)<<" , "<<c[j][y]<<" ]"<<endl;
		}
	}
	//Cmax=c[m-1][ind[n-1]];
	x = getElement(ind, n-1);
	Cmax = c[m-1][x];
	//cout<<"Cmax= "<<Cmax<<endl;
	return Cmax;
}
int NEH (int n, int m, list <list<int> > p, list<int> sigma1, list<int> sigma2, list<int> &S){
	//array des temps de traiement globeaux 
	int n1=sigma1.size();
	int n2=sigma2.size();
	int P[n2]={0};	int ind[n2];
	//meilleur Cmax  
	int C,C1,best;

	int a;
	// calcul des durees globales
	for(int i=0;i<n2;i++){
		//ind[i]=sigma2[i];
		a = getElement(sigma2, i);
		ind[i] = a;
		for(int j=0;j<m;j++) 	P[i]=P[i] + getElement(p, j, a);
	}
	//tri des taches selon l'ordre decroissant des duree globales
	Tri(n2,P,ind);// Tri est sans l'instruction ind[i]=i;
	// solution partielle
	list<int> s;
	S=sigma1;
	S.push_back(ind[0]);
	if(n2==1){
		C=ordonnancer(n,n,m,p,S);
	}else{
		for(int i=1; i<n2; i++){
		C1=n*m*5000;
		for(int k=0; k<=i; k++){
			s.clear(); 		s=S;
			list<int>::iterator v = s.begin();
			for(int w = 0; w<(n1+k); w++) v++;
			s.insert(v,ind[i]);
			//cout<<" "<<endl;
			//cout<<s.size()<<endl;
			C=ordonnancer(s.size(),n,m,p,s);
			if(C<C1){
				best=k;				C1=C;
			}	
		}
		S.insert(S.begin(),ind[i]);
	   }
	}
	
  //cout<<"S:"<<endl; for(int k=0; k<S.size();k++) cout <<"  "<<S[k]; cout<<" "<<endl;
	return C;
}


int LB (int n, int m, list <list<int> > p ){
	int p1[2][n];
	int p2[2][n];
	int L[n],r[n],q[n];
	int lb, minRq;  
	int c[2][n];
	int ind[n];
	int a;
	//construire l'instance de F2|li|Cmax
	lb=0;
	for(int k=0; k<m-1; k++){
		for(int l=k+1; l<m;l++){
		//	for(int i=0;i<n;i++){r[i]=0;	q[i]=0;		L[i]=0;	c[0][i]=0; c[1][i]=0; 	p1[0][i]=0; p1[1][i]=0;  p2[0][i]=0; p2[1][i]=0; }	
			for(int i=0; i<n; i++){
				//p1[0][i]=p[k][i];
				p1[0][i] = getElement(p, k, i);
				//p1[1][i]=p[l][i];
				p1[1][i] = getElement(p, l, i);
				L[i]=0; 
				if(l>k+1){
					for(int j=k+1; j<l; j++){
						//L[i]=L[i]+p[j][i];
						L[i] = L[i] + getElement(p, j, i);
					}
				}
				r[i]=0;
				if(k>0){
					for(int j=0; j<k; j++){
						//r[i]=r[i]+p[j][i];
						r[i]=r[i]+getElement(p,j,i);
					}
				}
				q[i]=0;
				if(l+1<m){
					for(int j=l+1; j<m; j++){
						//q[i]=q[i]+p[j][i];
						q[i]=q[i]+getElement(p,j,i);
					}	
				}
				//debut de resolution
				p2[0][i]=p1[0][i]+L[i];
				p2[1][i]=p1[1][i]+L[i];	
				ind[i]=i; // a ajouter 
			}//fin for (i)
			//resoudre le sous-probleme 
			for(int i=0;i<n;i++){
				for(int i1=i; i1<n; i1++){
					if(min(p2[0][i],p2[1][i1])>=min(p2[1][i],p2[0][i1])){
						permutation(p2[0][i],p2[0][i1]);
						permutation(p2[1][i],p2[1][i1]);
						permutation(ind[i],ind[i1]);
					}
				}
			}
			//calcul des temps de fin de traitement 
			c[0][ind[0]]=0+p1[0][ind[0]];
			c[1][ind[0]]=max(0,c[0][ind[0]])+L[ind[0]]+p1[1][ind[0]];
			for(int i=1; i<n; i++){
				c[0][ind[i]]=c[0][ind[i-1]]+p1[0][ind[i]];
				c[1][ind[i]]=max(c[1][ind[i-1]],c[0][ind[i]]+L[ind[i]])+p1[1][ind[i]];
			}
			//calcul du minRq
			minRq=n*m*10000;
			for(int i=0;i<n;i++){
				for(int i1=0;i1<n;i1++){
					if(i!=i1){
						if(r[i]+q[i1]<minRq){
							minRq=r[i]+q[i1];
						}
					}
				}
			}
			//cout<<"LB pour M"<<k<<" et M"<<l<<" est :"<< 
			//	c[1][ind[n-1]]<<" + "<<minRq<<" = "<<c[1][ind[n-1]]+minRq<<endl;
			if(lb<c[1][ind[n-1]]+minRq){
				lb=c[1][ind[n-1]]+minRq;
			}
		}//end for (l)
	}//end for (k)
	//cout<<" LB finale : "<<lb<<endl;
	return lb;
}//fin LB

int LB (int n, int m, list <list<int> >p, list<int> sigma1, list<int> sigma2 ){
	int n1=sigma1.size();
	int n2=sigma2.size();
	int ind[n2];
	int a,b,x;
	for(int i=0;i<n2;i++)	{
		ind[i]=getElement(sigma2,i);
	}	
	//ordonnancement de sigma1 
	int c[m][n];
	//initialiser la matrice "c" a 0
	for (int j=0;j<m;j++){	for(int i=0;i<n;i++) c[j][i]=0;		}
	//premiere tache
	for (int j=0;j<m;j++){
		for(int k=0;k<=j;k++){
			a = getElement(sigma1, 0);
			c[j][a]=c[j][a]+ getElement(p, k, a);
		}
	}
	//premiere machine 
	c[0][a]=0;
	for (int i=0;i<n1;i++){
		for(int k=0;k<=i;k++){
			a = getElement(sigma1, i);
			b = getElement(sigma1, k);
			c[0][a]=c[0][a]+getElement(p, 0, b);
		}
	}
	//le reste des taches et machines
	for(int i=1; i<n1;i++){
		for(int j=1; j<m;j++){
		    a =  getElement(sigma1, i);
		    b = getElement(sigma1, i-1);
			c[j][a]=max(c[j][b], c[j-1][a])+getElement(p, i, a);
			a = getElement(sigma1, i);
			b = getElement(sigma1, i-1);
			c[j][a] = max( c[j][b], c[j-1][a]) + getElement(p, j, a);
		}
	}
	//fin de l'ordonnancement de sigma1 
	//construction de F2|prmt, li|Cmax avec les taches de sigma2
	int p1[2][n2];
	int p2[2][n2];
	int L[n2],q[n2];
	int lb, minRq;  
	int c2[2][n2];
	lb=0;
	for(int k=0; k<m-1; k++){
		for(int l=k+1; l<m;l++){
			for(int i=0; i<n2; i++){
				a = getElement(sigma2, i);
				p1[0][i]=getElement(p,k,a);
				p1[1][i]=getElement(p,l,a);
				L[i]=0; 
				if(l>k+1){
					for(int j=k+1; j<l; j++)
						L[i]=L[i]+getElement(p, j, a);
				}//endif
				q[i]=0;
				if(l+1<m){
					for(int j=l+1; j<m; j++)
						q[i]=q[i]+getElement(p, j, a);
				}//endif
				//debut de resolution
				p2[0][i]=p1[0][i]+L[i];
				p2[1][i]=p1[1][i]+L[i];	
				ind[i]=i;
			}//fin for (i)
			//resoudre le sous-probleme 
			for(int i=0;i<n2;i++){
				for(int i1=i+1; i1<n2; i1++){
					if(min(p2[0][i],p2[1][i1])>=min(p2[1][i],p2[0][i1])){
						permutation(p2[0][i],p2[0][i1]);
						permutation(p2[1][i],p2[1][i1]);
						permutation(ind[i],ind[i1]);
					}
				}
			}
			//calcul des temps de fin de traitement du sous-probleme
			//ordonnancement de la premiere tache
			x = getElement(sigma1, n1-1);
			c2[0][ind[0]]=c[k][x]+p1[0][ind[0]];
			c2[1][ind[0]]=max(c2[0][ind[0]]+L[ind[0]],c[l][x])+p1[1][ind[0]];
			//ordonnancement du reste des taches
			for(int i=1; i<n2; i++){
				c2[0][ind[i]]=c2[0][ind[i-1]]+p1[0][ind[i]];
				c2[1][ind[i]]=max(c2[1][ind[i-1]],c2[0][ind[i]]+L[ind[i]])+p1[1][ind[i]];	
			}
			//Calcul de LB
			minRq=n*m*10000;
			for(int i=0;i<n2;i++){
				if(q[i]<minRq)
					minRq=q[i];
			}
			//cout<<"LB pour M"<<k<<" et M"<<l<<" est :"<< c2[1][ind[n2-1]]<<" + "<<
		 //minRq<<" = "<<c2[1][ind[n2-1]]+minRq<<endl;
			if(lb<c2[1][ind[n2-1]]+minRq){
				lb=c2[1][ind[n2-1]]+minRq;
			}
		}//end for (l)
	}//end for (k)
	//cout<<" LB finale : "<<lb<<endl;
	return lb;
}//fin LB(sigma1)

void GAP(int UB, int L, float& gap) {
	gap = float(UB - L);
	gap = gap / L;
	gap = gap * 100;
	if (gap == 0) {
		cout << it << "   " << L << "   " << UB << "   " << "0.00%" << endl;
	}
	else if (gap > 1) {
		cout << it << "   " << L << "   " << UB << "   "
			<< setprecision(3) << gap << "%" << endl;
	}
	else {
		cout << it << "   " << L << "   " << UB << "   "
			<< setprecision(2) << gap << "%" << endl;
	}

}

void print(list <int>& H) {
	int l = H.size();
	for (int i = 0; i < l; i++) {
		cout << " " << getElement(H, i) << "";
	}
	cout << " " << endl;
}

void Opt(int n, int m,  list <list <int> > p, bool& optimum, list<int> bestsol) {
	cout << "Optimum atient !!!!!" << endl;
	cout << "-----------------------------------------------------------"
		<< endl;
	cout << "La sequence d'ordrdonnancement optimale est : " << endl;
	print(bestsol);
	cout << "-----------------------------------------------------------"
		<< endl;
	ordonnancer(n,m,p,bestsol);
	optimum = true;
}

void BB(int n, int m, list <list <int> > p, list <int>& sigma,	list<int>& sigma1, 
	list <int>& sigma2, int& L, int& UB, list<int>& bestsol, bool& optimum, bool& stop) {
	if (stop == false) {
		/*
		cout<<"*************************"<<endl;
		
			cout<<"sigma : "<<endl;
			print(sigma);
			cout<<"sigma1 : "<<endl;
			print(sigma1);
			cout<<"sigma2 : "<<endl;
			print(sigma2);
			
		cout<<"*************************"<<endl;	
	*/
		float gap;
		it = it + 1;
		int j;
		list <int> sigma1_2;
		int UBsig1;
		
		int LBsigma1;
		if(sigma1.empty()==true){
			LBsigma1=L;
		}else{
			LBsigma1 = LB(n, m, p, sigma1, sigma2);
		} 
		//cout<<" LB sigma1 : "<<LBsigma1<<endl;
		//cout<<"LB : "<<L<<endl;
		if (LBsigma1 > L) {
			if (sigma.back() <= sigma2.size()) {
				j = sigma1.back();
				sigma2.insert(sigma2.begin(),j );				
				sigma1.pop_back();
			}
			else {
				while (sigma.back() == sigma2.size() + 1) {
					if (sigma.size() > 0) {
						if (sigma.size() == 1) {
							stop = true;//l'arbre a ete entierement parcouru 
						}
						j = sigma1.back();
						sigma2.insert(sigma2.begin(),j );						
						sigma1.pop_back();
						sigma.pop_back();
					}
				}
				if (sigma1.empty() == false) {
					j = sigma1.back();
					sigma2.insert(sigma2.begin(), j);					
					sigma1.pop_back();
				}
			}
		}
		else {
			UBsig1 = NEH(n,m, p, sigma1, sigma2, sigma1_2);
			if (UBsig1 == L) {
				UB = UBsig1;
				/**************** GAP ***************************************/
				GAP(UB, L, gap);
				bestsol.clear();
				bestsol = sigma1_2;
				Opt(n, m, p, optimum, bestsol);
				stop = true;
			}
			else {
				if (UBsig1 < UB) {
					UB = UBsig1;
					bestsol.clear();
					bestsol = sigma1_2;
					cout<<"la nouvelle solution : "<<endl;
					print(bestsol);
					/**************** GAP ***************************************/
					GAP(UB, L, gap);
					/**************** GAP ***************************************/
					cout << "sigma : ";
					print(sigma);
				}
			}
		}
		while (stop == false) {
			if (sigma1.size() == sigma.size()) {
				sigma.push_back(1);
			}
			else {
				j = sigma.back();
				j = j + 1;
				sigma.pop_back();
				sigma.push_back(j);
			}
			sigma1.push_back(sigma2.back());
			sigma2.pop_back();
			BB(n, m, p, sigma, sigma1, sigma2, L, UB, bestsol, optimum, stop);
		}
	}
}

int main(void){
	/*************************** INSTANCE *************************/
	cout<<"le nombre de taches ?"<<endl;
	int n;
	cin>>n;
	cout<<"le nombre de machines ?"<<endl;
	int m;
	cin>>m;
	list <list<int> > p;
	saisir(n, m, p);
	list <int> sigma; 
	list <int> sigma1;
	list <int> sigma2;
	for(int i=0; i<n; i++)
		sigma2.push_back(i);
	list <int> bestsol;
	bool optimum=false; 
	bool stop=false; 
	int UB=NEH(n,m,p,sigma1,sigma2, bestsol);
	cout<<" UB = "<<UB<<endl;
	int L=LB(n,m,p); 
	cout<<" LB = "<<L<<endl;
	sigma2.clear();
	
	sigma2=bestsol;
	cout<<"La sequence de branchement est :"<<endl;
	print(sigma2);
	float gap;
	if(L==UB){
		GAP(UB,L,gap);
		Opt(n,m,p,optimum, bestsol);
	}
	while(optimum==false){
		BB(n,m,p,sigma,sigma1,sigma2,L,UB,bestsol, optimum, stop);
		if(optimum==false){
			L=L+1;
			GAP(UB,L,gap);
			if(L==UB){
				Opt(n,m,p,optimum, bestsol);
			}else{
				stop=false;
			}
		}
	}

	
/*
2 4 1 6 3 5 8 
7 5 9 3 9 5 2 
5 6 3 9 7 1 6 
4 7 2 3 1 4 8

2 4 1 6 3 5 8 6 7 10 5 7 11 
7 5 9 3 9 5 2 4 8 9 12 4 9
5 6 3 9 7 1 6 8 2 6 9 4 10
4 7 12 3 2 6 8 2 11 4 9 7 4
*/
return 0;

}
	
	
	


	
