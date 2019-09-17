#include "euDataRead.h"

euDataRead::euDataRead(){
  
}

euDataRead::~euDataRead(){
}

bool euDataRead::open(string file){
  if (file.empty()){
    cout<< "euDataReader: Empty file name input." << endl; 
    return false;
  }
  string type_in = "raw"; // TODO: change to depend on the file name, more stable, MQ
  
  m_reader = eudaq::Factory<eudaq::FileReader>::MakeUnique(eudaq::str2hash(type_in), file);
  return true;
}

bool euDataRead::next (KpixEvent *kpixevent){
  auto ev = m_reader ->GetNextEvent();
  if (!ev){
    cout<< "euDataRead: empty event." << endl;
    return false;
  }
  
  std::vector<uint8_t> block = ev->GetBlock(0);
  size_t size_of_kpix = block.size()/sizeof(uint32_t);

  uint32_t *kpixdata = nullptr;
  if (size_of_kpix) // TODO: MQ, weird, I may have mistake with old code;
    kpixdata = reinterpret_cast<uint32_t *>(block.data());

  else {
    cout<< " " << endl;
    return false;
  }

  kpixevent->copy(kpixdata, size_of_kpix);
	 
  return true;
}

