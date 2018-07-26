  void plot(string descriptor, VD &field)
  {
    string filename = descriptor + ".plt";
    std::fstream f;
    f.open(filename.c_str(),std::ios::out);
    iLOOP 
      {
	jLOOP f << dx*(i-1) << " " << dx*(j-1) << " " << field[pid(i,j)] << endl;
	f << "\n";
      }
    f.close();
  }

  void plot(string descriptor, VD &xval, VD &yval, int len)
  {
    string filename = descriptor + ".plt";
    std::fstream f;
    f.open(filename.c_str(),std::ios::out);
    for ( int i = 1 ; i <= len ; ++i ) f << xval[i] << " " << yval[i] << endl;
    f.close();
  }
