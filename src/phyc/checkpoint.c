#include "checkpoint.h"
#include "filereader.h"

void checkpoint_apply(const char* file_path, Parameters* parameters) {
    StringBuffer* buffer = new_StringBuffer(10);
    Hashtable* hash = new_Hashtable_string(100);
    hashtable_set_key_ownership(hash, true);
    hashtable_set_value_ownership(hash, true);
    FileReader* reader = new_FileReader(file_path, 1000);
    while (reader->read_line(reader)) {
        StringBuffer_trim(reader->buffer);
        if (reader->buffer->length == 0) {
            continue;
        }

        int c = 0;
        const char* pch = reader->line;
        while (*pch != '\0') {
            if (*pch == ',') {
                break;
            }
            c++;
            pch++;
        }
        StringBuffer_set_nstring(buffer, reader->line, c);
        Hashtable_add(hash, String_clone(buffer->c), new_Double(atof(pch + 1)));
    }
    free_StringBuffer(buffer);
    free_FileReader(reader);

    for (size_t i = 0; i < Parameters_count(parameters); i++) {
        double* value = Hashtable_get(hash, Parameters_name(parameters, i));
        if (value != NULL) {
            Parameters_set_value(parameters, i, *value);
        }
    }
    free_Hashtable(hash);
}

void checkpoint_save(const char* file_path, Parameters* parameters) {
    bool overwrite = file_exists(file_path);
    StringBuffer* buffer = new_StringBuffer(10);

    if (overwrite) {
        StringBuffer_set_string(buffer, file_path);
        do {
            StringBuffer_append_string(buffer, ".new");
        } while (file_exists(buffer->c));
        file_copy(file_path, buffer->c);
    }

    FILE* file = fopen(file_path, "w");
    for (size_t i = 0; i < Parameters_count(parameters); i++) {
        Parameter* p = Parameters_at(parameters, i);
        fprintf(file, "%s,%e\n", Parameter_name(p), Parameter_value(p));
    }
    fclose(file);

    if (overwrite) {
        remove(buffer->c);
        free_StringBuffer(buffer);
    }
}