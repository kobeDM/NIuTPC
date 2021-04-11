// -*- C++ -*-

#ifndef DCEVENT_INCLUDED_H
#define DCEVENT_INCLUDED_H

#include <assert.h>
#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include <netdb.h>
#include <queue>

//
class dc_header {
public:
    unsigned int m_magic;
    unsigned int m_id;
    unsigned int m_length;
    unsigned int m_trig_cnt;

    unsigned int magic() const {
        return ntohl(m_magic);
    };

    unsigned int id() const {
        return ntohl(m_id);
    };

    unsigned int length() const {
        return ntohl(m_length);
    };

    unsigned int trig_cnt() const {
        return ntohl(m_trig_cnt);
    };

    void show() const {
        printf("magic = %08x, id = %08x, length = %08x, trig cnt = %08x\n", 
	       magic(),id(),length(),trig_cnt());
    };
};

//
class dc_event {
private:
    dc_header m_header;
    unsigned char* m_body;
    int m_len;

public:
    dc_event(dc_header& header, unsigned char* p, int len) 
	: m_header(header), m_len(len)  {
	
	m_body = new unsigned char [len];
	assert(m_body);
        //printf("p = %p, len = %d\n", p, m_len);

        memcpy(m_body, p, len);
    }

    dc_event(const dc_event & e) {
	m_header = e.m_header;
        m_len    = e.m_len;
        
        m_body = new unsigned char[m_len];

        memcpy(m_body, e.m_body, m_len);
    };

    ~dc_event() {
        delete [] m_body;
    };


    void show() const {
        printf("=============dc_event.show()=======================\n");

        m_header.show();
        assert(m_len % 4 == 0);

        int wlen = m_len / 4;
        unsigned int * p = (unsigned int *)m_body;

        for (int i=0; i<wlen; i++) {
            if (i % 4 == 0)
                printf("%04x", i*4);
            printf(" %08x", ntohl(p[i]));
            if (i % 4 == 3)
                fputc('\n', stdout);
        }

        if (wlen % 4 != 0)
            fputc('\n', stdout);

        printf("====================================\n");
    };

    const dc_header & header() const {
        return m_header;
    };

    int len() const {
        return m_len;
    };

    unsigned char * body() const {
        return m_body;
    };
};


#endif  // DCEVENT_INCLUDED_H
