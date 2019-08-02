//-----------------------------------------------------------------------------
// File          : ymlReader.cxx
// Author        : Mengqing Wu <mengqing.wu@desy.de>
// Created       : 05/10/2018
// Project       : Similar to the old XmlVariables, a new Yaml parser/handler.
//-----------------------------------------------------------------------------

#include <iostream>
//#include <iomanip>
#include <stdio.h>
//#include <stdarg.h>
#include <math.h>
#include <fstream>
#include <string.h>

#include "TFile.h"

#include <yaml.h>
#include <assert.h>

using namespace std;


void compactParser(const char* fn);

int exParser1(const char* fn);
int exParser2(const char* fn);

int main(int argc, const char* argv[]){

  if (argc<2){
    cout << "Usage: ymlReader config.yml"<<endl;
    return 0;
  }
 
  compactParser(argv[1]);
  //  exParser1(argv[1]);
  
  return 1;
}

// Process yml
void compactParser(const char* fn){

  /* funcs:
   * 1) a fake DataRead wrapper;
   * 2) 
   */
  
  std::ifstream ifs;
  
  ifs.open(fn, std::ifstream::in);
  if (!ifs.is_open()) {

    cout<< "Error: in file not open" << endl;
    return;
  }
  
  //  const auto offset = 100;
  //  ifs.seekg(offset);
  int length = 100;
  char* pChars = new char[100];

  ifs.seekg(0, ios::beg);
  
  // read the data chunk to buffer and from buffer to output stream
  //std::vector buff (100, 0);
  //ifs.read(buff.data(), buff.size());
  //os.write(buff.data(), buff.size());

  ifs.read(pChars, length);


  // -- debug START: output file to test
  std::ofstream os;
  os.open("test.txt", ios::out);
  os.write(pChars, length);
  os.close();
  // -- debug END
  
  ifs.close();
}


int exParser2(const char* fn)
{

  FILE *file = fopen(fn, "r");
  yaml_parser_t parser;
  yaml_document_t document;
  yaml_node_t *node;
  int i = 1;

  assert(file);
  assert(yaml_parser_initialize(&parser));

  // you can use yaml_parser_set_input_string for string input
  yaml_parser_set_input_file(&parser, file);
  
  if (!yaml_parser_load(&parser, &document)) {
    goto done;
  }

  while(1) {
    node = yaml_document_get_node(&document, i);
    if(!node) break;
    printf("Node [%d]: %d\n", i++, node->type);
    if(node->type == YAML_SCALAR_NODE) {
      printf("Scalar [%d]: %s\n", node->data.scalar.style, node->data.scalar.value);
    }
  }
  yaml_document_delete(&document);


 done:
  yaml_parser_delete(&parser);
  assert(!fclose(file));
  
  return 1;
}


int exParser1(const char* fn){

  FILE *fh = fopen( fn, "r");
  yaml_parser_t parser;
  yaml_token_t  token;   /* new variable */
  //  yaml_event_t  event;   /* New variable */
    
  /* Initialize parser */
  if(!yaml_parser_initialize(&parser))
    fputs("Failed to initialize parser!\n", stderr);
  if(fh == NULL)
    fputs("Failed to open file!\n", stderr);

  /* Set input file */
  yaml_parser_set_input_file(&parser, fh);

  /* CODE HERE */

  /* BEGIN new code */
  do {
    yaml_parser_scan(&parser, &token);
    switch(token.type)
    {
    /* Stream start/end */
    case YAML_STREAM_START_TOKEN: puts("STREAM START"); break;
    case YAML_STREAM_END_TOKEN:   puts("STREAM END");   break;
    /* Token types (read before actual token) */
    case YAML_KEY_TOKEN:   printf("(Key token)   "); break;
    case YAML_VALUE_TOKEN: printf("(Value token) "); break;
    /* Block delimeters */
    case YAML_BLOCK_SEQUENCE_START_TOKEN: puts("<b>Start Block (Sequence)</b>"); break;
    case YAML_BLOCK_ENTRY_TOKEN:          puts("<b>Start Block (Entry)</b>");    break;
    case YAML_BLOCK_END_TOKEN:            puts("<b>End block</b>");              break;
    /* Data */
    case YAML_BLOCK_MAPPING_START_TOKEN:  puts("[Block mapping]");            break;
    case YAML_SCALAR_TOKEN:  printf("scalar %s \n", token.data.scalar.value); break;
    /* Others */
    default:
      printf("Got token of type %d\n", token.type);
    }
    if(token.type != YAML_STREAM_END_TOKEN)
      yaml_token_delete(&token);
  } while(token.type != YAML_STREAM_END_TOKEN);
  yaml_token_delete(&token);
  /* END new code */
  

  
  /* Cleanup */
  yaml_parser_delete(&parser);
  fclose(fh);
  return 0;
}
