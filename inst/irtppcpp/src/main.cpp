#include <iostream>
#include <estimation/emestimation.h>
#include <utils/Input.h>

using namespace irtpp;
using std::cout;
using std::endl;

int main(int argc, char ** argv)
{
  void**   status_list;
  dataset* d;
  Input input;

  d = new dataset(0);
  input.importCSV(argv[1], *d, 1, 0);

  emestimation em1(new onepl(),   d);
  //emestimation em2(new twopl(),   d);
  //emestimation em3(new threepl(), d);

  status_list = em1.estimate();
  cout << "Iterations: " << *((int*)status_list[0]) << endl;
  delete (int*)status_list[0];
  delete (bool*)status_list[1];
  delete [] status_list;
  // cout << "------------" << endl;
  // status_list = em2.estimate();
  // cout << "Iterations: " << *((int*)status_list[0]) << endl;
  // delete (int*)status_list[0];
  // delete (bool*)status_list[1];
  // delete [] status_list;
  // cout << "------------" << endl;
  // status_list = em3.estimate();
  // cout << "Iterations: " << *((int*)status_list[0]) << endl;
  // delete (int*)status_list[0];
  // delete (bool*)status_list[1];
  // delete [] status_list;

  delete d;

  return 0;
}