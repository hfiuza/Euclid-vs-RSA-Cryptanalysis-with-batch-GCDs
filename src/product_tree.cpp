#include <iostream>
#include <vector>
using namespace std;

long long product(vector<long long>set_of_numbers){
    int n = set_of_numbers.size();
    if(n==0)
        return 1;
    else if(n==1)
        return set_of_numbers[0];
    else{
        for(int i=0; i<n; i+=2){
            if(i!=(n-1))
                set_of_numbers[i/2] = set_of_numbers[i]*set_of_numbers[i+1]; 
            else
                set_of_numbers[i/2] = set_of_numbers[i];
        }
        set_of_numbers.resize((n+1)/2);    
        return product(set_of_numbers);
    }   
}

vector<vector<long long> > product_tree(vector<long long>set_of_numbers){
    vector<vector<long long> >result;
    result.push_back(set_of_numbers);
    int n = set_of_numbers.size();
    while(n>1){
         for(int i=0; i<n; i+=2){
            if(i!=(n-1))
                set_of_numbers[i/2] = set_of_numbers[i]*set_of_numbers[i+1]; 
            else
                set_of_numbers[i/2] = set_of_numbers[i];
        }
        n = (n+1)/2;
        set_of_numbers.resize(n);    
        result.push_back(set_of_numbers); 
    }
    return result;
}
void print_product_tree(vector<vector<long long> >p_T){
    int k  =p_T.size();
    for(int i=0; i<k; i++){
        cout << "[ ";
        for(int j=0; j<p_T[i].size(); j++){
            cout << p_T[i][j] << " ";
        }
        cout << "]" << endl;
    }
    cout << endl;
}
vector<long long> remainders_using_product_tree(long long M, vector<vector<long long> >p_T){
    vector<long long> result(1, M);
    int k = p_T.size();
    int len = 1; // represents the current size of result = size of p_T[i]
    for(int i=k-1; i>=0; i--){
        len =  p_T[i].size();
        result.resize(len);
        for(int j=len-1; j>=0; j--){
            result[j] = result[j/2] % p_T[i][j];   
        }
    }
    return result;
}
vector<long long> remainders(long long M, vector<long long>set_of_numbers){
    vector<vector<long long> > my_tree = product_tree(set_of_numbers);
    return remainders_using_product_tree(M, my_tree);
}
void print_remainders(vector<long long> rM){
    int k = rM.size();
    cout << "[ ";
    for(int i=0; i<k; i++){
        cout << rM[i] << " ";
    }
    cout << "]" << endl;
}
vector<long long>square_numbers(vector<long long>&set_of_numbers){
    int k = set_of_numbers.size();
    vector<long long>sq_n;
    for(int i=0; i<k; i++){
        sq_n.push_back(set_of_numbers[i]*set_of_numbers[i]);
    }
    return sq_n;
}

long long gcd(long long a, long long b){
    if(a==0)
        return b;
    else
        return gcd(b%a, a);
}

vector<long long>compute_gcds(vector<long long>&set_of_numbers){
    vector<long long>squared_numbers = square_numbers(set_of_numbers);
    vector<long long>rems = remainders(product(set_of_numbers), squared_numbers);
    vector<long long>gcds;
    int k = set_of_numbers.size();
    for(int i=0; i<k; i++)
        gcds.push_back(gcd(set_of_numbers[i], rems[i]/set_of_numbers[i]));        
    return gcds;
}
vector<long long>gcds_treatment(vector<long long>&set_of_numbers, vector<long long>&gcds){
    vector<long long>factorizations;
    int k = set_of_numbers.size();
    for(int i=0; i<k; i++){
        if(gcds[i]==1) // we can't break the key
            factorizations.push_back(0);
        if(gcds[i] < set_of_numbers[i]){
            if (gcds[i]*gcds[i] <= set_of_numbers[i])
                factorizations.push_back(gcds[i]);
            else
                factorizations.push_back(set_of_numbers[i]/gcds[i]);
        }
        else{
            bool is_repeated=false;
            for(int j=0; j<k; j++){
                if(set_of_numbers[j]==set_of_numbers[i] && i!=j){
                    is_repeated = true;
                    break;
                }
            }
            if(is_repeated)
                factorizations.push_back(0);
            else{
                bool pushed = false;
                for(int j=0; j<k; j++){
                    long long temp = gcd(set_of_numbers[i], set_of_numbers[j]);
                    if(j!=i && temp>1){
                        if(temp*temp <= set_of_numbers[i])
                            factorizations.push_back(temp);
                        else
                            factorizations.push_back(set_of_numbers[i]/temp);
                        pushed = true;
                        break;
                    }
                }
                if(!pushed)
                    cout << "ERROR ERROR ERROR " << endl;
            }
            
        }
    }
    return factorizations;
}


int main(){
    long long my_nbs[] = {21, 35, 707, 221, 33, 1843};
    vector<long long>my_numbers(my_nbs, my_nbs+sizeof(my_nbs)/sizeof(long long));
    //cout << product(my_numbers) << endl;
    //print_product_tree(product_tree(my_numbers));
    //print_remainders(remainders(1001, my_numbers));   
    
    vector<long long>gcds = compute_gcds(my_numbers);
    vector<long long>factors = gcds_treatment(my_numbers, gcds);
    for(int i=0; i<factors.size(); i++){
        cout << factors[i] << " ";
    }
    cout << endl;
    
        
}



