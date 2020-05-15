#include <iostream>

using namespace std;

int main(){
  // cout << "hello world!" << endl;
  //
  //
  // int x;
  //
  // cout << "type int:";
  // cin >> x;
  //
  // x *= 2;
  //
  // cout << x << endl;
  //
  // string s;
  // cout << "type string:";
  // cin >> s;
  //
  // int len = s.size();
  // string s2;
  // for(int i=0; i <= len; i++){
  //   s2 += s[len-i];
  //
  // }
  //
  // cout << s2 << endl;
  //
  //
  // int y = 5;
  //
  // cout << ++y << " " << y << endl;
  //
  // cout << y++ << " " << y << endl;


  int a = 0;
  int b[7564];
  try {
    cout << b[a++] << endl;
  }
  catch (){
    cout << "lenght of array is " << a << endl;
  }


  return 0;
}
