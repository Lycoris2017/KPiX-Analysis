//-----------------------------------------------------------------------------
// File          : ymlReader.cxx
// Author        : Mengqing Wu <mengqing.wu@desy.de>
// Created       : 19/10/2018
// Project       : Write / emitt what you read from memory
//-----------------------------------------------------------------------------

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <string.h>
#include <map>

#include <yaml.h>

//#include <stdlib.h>

// Define variable holder
typedef std::map<std::string,std::string> varHolder; // <name, value>

varHolder data = {
  {"CalState", "Inject"},
  {"CalChannel", "0"},
  {"CalDac", "200"},
};


int main(int argc, char *argv[])
{
    yaml_emitter_t emitter;
    yaml_event_t event;

    struct fruit *f;
    char buffer[64];

    yaml_emitter_initialize(&emitter);
    yaml_emitter_set_output_file(&emitter, stdout);

    //FILE *output = fopen("out_emitter.yml","wb");
    //yaml_emitter_set_output_file(&emitter, output);
 
    yaml_stream_start_event_initialize(&event, YAML_UTF8_ENCODING);
    if (!yaml_emitter_emit(&emitter, &event)) goto error;
 
    yaml_document_start_event_initialize(&event, NULL, NULL, NULL, 0);
    if (!yaml_emitter_emit(&emitter, &event)) goto error;
 
    yaml_sequence_start_event_initialize(&event, NULL, (yaml_char_t *)YAML_SEQ_TAG,
					 1, YAML_ANY_SEQUENCE_STYLE);
    if (!yaml_emitter_emit(&emitter, &event)) goto error;

    //-- start to fill the map:
    for (auto & kv : data) {
      //std::cout << kv.first << ": " << kv.second << std::endl;

      //-- start map tag
      yaml_mapping_start_event_initialize(&event, NULL, (yaml_char_t *)YAML_MAP_TAG,
					  1, YAML_ANY_MAPPING_STYLE);
      if (!yaml_emitter_emit(&emitter, &event)) goto error;

      //-- fill token key
     
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
				   (yaml_char_t *)kv.first.c_str(), kv.first.size(),
				   1, 0, YAML_PLAIN_SCALAR_STYLE);
      if (!yaml_emitter_emit(&emitter, &event)) goto error;
 
      //-- fill token value
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
				   (yaml_char_t *)kv.second.c_str(), kv.second.size(),
				   1, 0, YAML_PLAIN_SCALAR_STYLE);
      if (!yaml_emitter_emit(&emitter, &event)) goto error;
            
      //-- end map tag
      yaml_mapping_end_event_initialize(&event);
      if (!yaml_emitter_emit(&emitter, &event)) goto error;
      
	
    }
 
    yaml_sequence_end_event_initialize(&event);
    if (!yaml_emitter_emit(&emitter, &event)) goto error;
 
    yaml_document_end_event_initialize(&event, 0);
    if (!yaml_emitter_emit(&emitter, &event)) goto error;
 
    yaml_stream_end_event_initialize(&event);
    if (!yaml_emitter_emit(&emitter, &event)) goto error;
 
    yaml_emitter_delete(&emitter);
    return 0;
 
error:
    fprintf(stderr, "Failed to emit event %d: %s\n", event.type, emitter.problem);
    yaml_emitter_delete(&emitter);
    return 1;
}
