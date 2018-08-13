  void plot(string descriptor, VD &field)
  {
    string filename = Name + "_" + descriptor + ".plt";
    std::fstream f;
    f.open(filename.c_str(),std::ios::out);
    iLOOP 
      {
	jLOOP f << dx*(i-1) << " " << dx*(j-1) << " " << field[pid(i,j)] << endl;
	f << "\n";
      }
    f.close();
  }

