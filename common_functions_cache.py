
# coding: utf-8

# In[1]:

'''Version of common functions with more caching, also using new sqlite db schema (LTCC_cache). Functions for querying abstracts on ePMC, saving results in sqlite db and displaying results in dataframe.
More checks in this version to avoid duplicate retrieval of information from web services. 
Functions available:
create_db(db_name)
pop_chembl_pmids(db_name) -- this populates the chembl_pmids table with pmids from a specific chembl_version.
def define_synonym_queries(term_dict_1, term_dict_2 = None, term_dict_3 = None)
get_hit_profile(query_list)
get_pmids(query, db_name) -- query is string
get_article_data(query_id, db_name) -- query is string
get_pmids_and_article_data(query, db_name) -- query is string
get_availabilities(query_id, db_name)
get_scores(query_id, db_name)
set_chembl_values(query_id, db_name)
get_df(query_id_list, db_name, sql_condition=None)
separate_column_df(query_id)
colour_terms(df, markup_list)
plot_scores(query_id_list, db_name)
'''


# In[ ]:

import sqlite3 as lite
import datetime
import requests
from lxml import etree
import cx_Oracle
import re
import shelve
import base64
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import sqlite3 as lite
import datetime
from IPython.display import HTML
import lxml.html
from time import sleep


# In[3]:

def create_db(db_name):
    '''Create the SQLite database in the current directory.
    kwargs: db_name -- name of the new SQLite db you are creating'''
    
    conn = lite.connect(db_name)
    cursor = conn.cursor()
    
    #create various tables
    cursor.execute("create table queries(query_id integer primary key, query text, hitcount integer, date_performed text)")
    cursor.execute("create table result_ids(query_id integer, pmid integer, primary key(query_id, pmid))")
    cursor.execute("create table article_data(pmid integer primary key, year integer, title text, abstract text, journal_title text, journal_abbrev_title text, in_epmc integer, avail_codes text, pdf_links text, other_links text, in_chembl int)")
    cursor.execute("create table article_links(pmid integer primary key, campus_links text, request_access text)")
    
    cursor.execute("create table scores(pmid integer primary key, score real)")
    cursor.execute("create table chembl_pmids(pmid integer primary key)")
    
    cursor.execute("create table error_records(query_id integer, object_id text, pmid integer, error_comment text)")
    
    conn.commit()
    conn.close()
    
    return None


# In[ ]:

def pop_chembl_pmids(db_name):
    '''Get all pubmed ids from ChEMBL_20 and populate the chembl_pmids table in the SQLite db with those pmids. Needs a file with login details to ChEMBL database.
    kwargs: db_name -- name of the SQLite db of which the chembl_pmids table should be updated'''
    
    # get login details from file
    fileObj = open('/homes/ines/chembl_20_login_details.txt')
    database = fileObj.read().strip()
    fileObj.close()
    
    conn = cx_Oracle.connect(database)
    cursor = conn.cursor()
    chembl_pmid_list = []
    
    # get all pmids in chembl_20. Need to redo this for new version of chembl
    sql = '''
    select distinct pubmed_id from docs where pubmed_id is not null
    '''
    cursor.execute(sql)
    
    chembl_pmids_list = [pmid_tuple[0] for pmid_tuple in cursor.fetchall()]

    conn.close()
    
    # open queries_db
    conn = lite.connect(db_name)
    cursor = conn.cursor()
    
    # insert pmids in chembl_pmids table
    for pmid in chembl_pmids_list:
        cursor.execute("insert or ignore into chembl_pmids(pmid) values(?)", (pmid,))
    
    conn.commit()
    conn.close()
    
    return None


# In[4]:

def define_synonym_queries(term_dict_1, term_dict_2 = None, term_dict_3 = None):
    """Submit one, two or three dicts with synonyms/alternative terms to be used. 
    If there are no synonyms bu just one term, just include the primary term as one item in the list that is that value of the dict.
    For each dict the terms in the dict will be joined by OR. The three different dicts will be joined by AND in the query.
    Returns a list of queries. It is possible to have only one item in a dictionary, if you only want to do one query.
    Otherwise include the independent items to look for as different items in the dictionary.
    
    kwargs: term_dict_1 -- dict of terms (in list in value of dict) to be connected with OR in the query
            term_dict_2 -- additional terms to be included in the same query (default = None)
            term_dict_3 -- even more terms to be included in the same query (default = None)
    """

    my_queries = []
    
    if term_dict_2 == None:
    
        for term in term_dict_1:
            query_terms = '"'+'" OR ABSTRACT:"'.join(term_dict_1[term])+'"'
            query = '(ABSTRACT:'+query_terms+')'
            my_queries.append(query)

    elif term_dict_2 != None and term_dict_3 == None:
        
        for term in term_dict_1:
            for term_2 in term_dict_2:
                    
                query_dict_1 = '"'+'" OR ABSTRACT:"'.join(term_dict_1[term])+'"'
                query_dict_2 = '"'+'" OR ABSTRACT:"'.join(term_dict_2[term_2])+'"'
                query = '(ABSTRACT:'+query_dict_1+')'+' AND '+'(ABSTRACT:'+query_dict_2+')'
                my_queries.append(query)
        
    elif term_dict_3 != None:
        
        for term in term_dict_1:
            for term_2 in term_dict_2:
                for term_3 in term_dict_3:
                    
                    query_dict_1 = '"'+'" OR ABSTRACT:"'.join(term_dict_1[term])+'"'
                    query_dict_2 = '"'+'" OR ABSTRACT:"'.join(term_dict_2[term_2])+'"'
                    query_dict_3 = '"'+'" OR ABSTRACT:"'.join(term_dict_3[term_3])+'"'
                    query = '(ABSTRACT:'+query_dict_1+')'+' AND '+'(ABSTRACT:'+query_dict_2+')'+' AND '+'(ABSTRACT:'+query_dict_3+')'
                    my_queries.append(query)
        
    return my_queries


# In[ ]:

def get_hit_profile(query_list):
    
    """For each query in query_list print profile hits. This uses the profile module of ePMC webservices. Will print results to console.
    kwargs: query_list -- list of queries (strings), e.g from define_synonym_queries function."""
    
    st_dict = {}
    base = 'http://www.ebi.ac.uk/europepmc/webservices/rest/profile/query=({} AND (src:MED OR src:PMC OR src:CTX))'
    total_nr_articles = 0
    
    for query in query_list:
            
        response = requests.get(base.format(query))
        tree = etree.fromstring(response.content)
        
        if response.status_code != 200:
            print('for query:{} status_code: {}'.format(query, str(response.status_code)))
            
        all_result = tree.xpath('/responseWrapper/profileList/pubType[@name="ALL"]')
        ft_result = tree.xpath('/responseWrapper/profileList/pubType[@name="FULL TEXT"]')
            
        all_string = str(etree.tostring(all_result[0], encoding = 'unicode'))
        ft_string = str(etree.tostring(ft_result[0], encoding = 'unicode'))
                
        if re.search(r'count="0"', all_string) == None:
            st_dict[query] = {'all':int(float(re.search(r'count="([0-9]+)"', all_string).group(1))), 'full text':re.search(r'count="([0-9]+)"', ft_string).group(1)}
        else:
            continue
    
    for entry in st_dict:
        print(entry)
        print(st_dict[entry])
        total_nr_articles = total_nr_articles + st_dict[entry]['all']
    
    print('\n'+'total number of hits = '+str(total_nr_articles))
    
    return None


# In[5]:

def get_pmids(query, db_name):
    '''For given query (string) get results in forms of idlist from ePMC. Record the query and date in the queries table. Save pmids associated with the query in the result_ids table.
    It will be better to use the get_pmids_and_article data function, but if only getting the pmids is required for quick overlap checking, for example, then use this function.
    Return the query_id assigned to the query in the queries table, can then be used in subsequent functions.
    kwargs: query -- should be string
            db_name
    '''
    
    moment = datetime.datetime.now()
    current_date = "{}-{}-{}".format(moment.year, moment.month, moment.day)
    
    base = 'http://www.ebi.ac.uk/europepmc/webservices/rest/search/query=({} AND (src:MED OR src:PMC OR src:CTX))&resulttype=idlist&page={}'
    
    conn = lite.connect(db_name)    
    cursor = conn.cursor()
    
    def do_query(my_query, page_nr = 1):
        
        response = requests.get(base.format(my_query, page_nr))
        
        assert response.status_code == 200
        
        tree = etree.fromstring(response.content)
        
        return tree
    
    def save_pmid(result):
        
        pmid = result.xpath('pmid/text()')[0]
        
        cursor.execute("insert or ignore into result_ids(query_id, pmid) values (?,?)", (current_query_id, pmid))
        
        return None
    
    
    tree = do_query(query)
 
    result_hitcount = int(tree.xpath('/responseWrapper/hitCount/text()')[0]) 
    cursor.execute("insert into queries(query_id, query, hitcount, date_performed) values (NULL,?,?,?)", (query, result_hitcount, current_date))
    current_query_id = cursor.lastrowid


    result_pmids = tree.xpath('/responseWrapper/resultList/result')


    for item in result_pmids:

        try:
            save_pmid(item)
        except IndexError:
            object_id = item.xpath('id/text()')[0]
            error_comment = '(get_pmids) - IndexError with XML, possibly no pmid for this item or no journal info e.g. when is book chapter'
            cursor.execute("insert or ignore into error_records(query_id, object_id, error_comment) values (?,?,?)", (current_query_id, object_id, error_comment))
            continue


    if (result_hitcount % 25 == 0) and (result_hitcount / 25 > 1): # hitCount is multiple of 25 (remainder==0) but more than one page
        total_pages = int(result_hitcount / 25)

    else:  
        total_pages = (int(result_hitcount // 25)) + 1

    for page in range(2,total_pages+1): #page one was already done so continue from page 2 and add results to same dictionary

        tree = do_query(query, page_nr = page)
        result_pmids = tree.xpath('/responseWrapper/resultList/result')

        for item in result_pmids:

            try:
                save_pmid(item)
            except IndexError:
                object_id = item.xpath('id/text()')[0]
                error_comment = '(get_pmids) - IndexError with XML, possibly no pmid for this item or no journal info e.g. when is book chapter'
                cursor.execute("insert or ignore into error_records(query_id, object_id) values (?,?,?)", (current_query_id, object_id, error_comment))
                continue
                        
    conn.commit()
    conn.close()

    return current_query_id


# In[4]:

def get_article_data(query_id, db_name):
    '''If already perfomed the get_pmids function for a query, now go in more detail and get the article data associated with the pmids form ePMC.
    Checks whether the pmid already exists in the article_data table and if so, skips that one and does not extract article data from ePMC.
    Returns None. 
    kwargs: query_id -- (obtained from get_pmids function)
            db_name'''
    
    conn = lite.connect(db_name)    
    cursor = conn.cursor()
    
    base = 'http://www.ebi.ac.uk/europepmc/webservices/rest/search/query=({} AND (src:MED OR src:PMC OR src:CTX))&resulttype=core&page={}'
    
    cursor.execute('select query from queries where query_id = ?', (query_id,))
    query = cursor.fetchall()[0][0]
    
    

    def do_query(my_query, page_nr = 1):

        response = requests.get(base.format(my_query, page_nr))
        assert response.status_code == 200, "status code != 200"
    
        tree = etree.fromstring(response.content)
        return tree

    def retrieve_data(result):
    
        #pmid = result.xpath('pmid/text()')[0]
        title = result.xpath('title/text()')[0]
        abstract = result.xpath('abstractText/text()')[0]
        year = result.xpath('journalInfo/yearOfPublication/text()')[0]
        journal_title = result.xpath('journalInfo/journal/title/text()')[0]
        journal_abbrev_title = result.xpath('journalInfo/journal/medlineAbbreviation/text()')[0]
        
        if result.xpath('inEPMC/text()')[0] == 'Y': 
            in_ePMC = 1 
        else: 
            in_ePMC = 0

        url_list = []
        avail_code_list = []
        doc_style_list = []
        other_links_list = []
        pdf_links_list = []

        ft_url_list = result.xpath('fullTextUrlList/fullTextUrl')
        
        for url_item in ft_url_list:
            avail_code = url_item.xpath('availabilityCode/text()')[0]
            avail_code_list.append(avail_code)
            doc_style = url_item.xpath('documentStyle/text()')[0]
            doc_style_list.append(doc_style)

            if avail_code != 'S' and doc_style != 'pdf':
                url = url_item.xpath('url/text()')[0]
                other_links_list.append(url)

        if 'pdf' in doc_style_list:
            for url_item in ft_url_list:
                doc_style = url_item.xpath('documentStyle/text()')[0]
                if doc_style == 'pdf':
                    url = url_item.xpath('url/text()')[0]
                    pdf_links_list.append(url)                        
        else:
            pass

        avail_codes = ', '.join(set(avail_code_list)) if avail_code_list else None
        other_links = ', '.join(other_links_list) if other_links_list else None
        pdf_links = ', '.join(pdf_links_list) if pdf_links_list else None
        
        cursor.execute("insert or ignore into article_data(pmid, year, title, abstract, journal_title, journal_abbrev_title, in_epmc, avail_codes, pdf_links, other_links) values (?,?,?,?,?,?,?,?,?,?)", (pmid, year, title, abstract, journal_title, journal_abbrev_title, in_ePMC, avail_codes, pdf_links, other_links))
        
        conn.commit()
        
        return None
    
    # now start the tasks
    
    cursor.execute('select r.pmid from result_ids r where query_id = ? and r.pmid in (select distinct a.pmid from article_data a)', (query_id,))
    existing_pmid_list = [i[0] for i in cursor.fetchall()]

    tree = do_query(query)
    
    result_abstracts = tree.xpath('/responseWrapper/resultList/result')
    result_hitcount = int(tree.xpath('/responseWrapper/hitCount/text()')[0]) 

    for item in result_abstracts:
        try:
            pmid = int(item.xpath('pmid/text()')[0])

            if pmid in existing_pmid_list:
                continue

            else:
                try:
                    retrieve_data(item)
                except IndexError:
                    object_id = item.xpath('id/text()')[0]
                    error_comment = '(get_article_data) - IndexError with XML, possibly field not present, e.g. no journal info e.g. when is book chapter'
                    cursor.execute("insert or ignore into error_records(query_id, object_id, error_comment) values (?,?,?)", (current_query_id, object_id, error_comment))
                    continue
        except IndexError:
            object_id = item.xpath('id/text()')[0]
            error_comment = '(get_article_data) - IndexError with XML, possibly no pmid for this item'
            cursor.execute("insert or ignore into error_records(query_id, object_id, error_comment) values (?,?,?)", (current_query_id, object_id, error_comment))
            continue
        
    if (result_hitcount % 25 == 0) and (result_hitcount / 25 > 1): # hitCount is multiple of 25 (remainder==0) but more than one page
        total_pages = int(result_hitcount / 25)

    else:  
        total_pages = (int(result_hitcount // 25)) + 1 # this holds for cases where there is only one page as well
        
    for page in range(2,total_pages+1): #page one was already done so continue from page 2 and add results to same dictionary

        tree = do_query(query, page_nr = page)
        result_abstracts = tree.xpath('/responseWrapper/resultList/result')

        for item in result_abstracts:
            try:
                pmid = int(item.xpath('pmid/text()')[0])

                if pmid in existing_pmid_list:
                    continue

                else:
                    try:
                        retrieve_data(item)
                    except IndexError:
                        object_id = item.xpath('id/text()')[0]
                        error_comment = '(get_article_data) - IndexError with XML, possibly field not present, e.g. no journal info e.g. when is book chapter'
                        cursor.execute("insert or ignore into error_records(query_id, object_id, error_comment) values (?,?,?)", (current_query_id, object_id, error_comment))
                        continue
            except IndexError:
                object_id = item.xpath('id/text()')[0]
                error_comment = '(get_article_data) - IndexError with XML, possibly no pmid for this item'
                cursor.execute("insert or ignore into error_records(query_id, object_id, error_comment) values (?,?,?)", (current_query_id, object_id, error_comment))
                continue


    conn.commit()
    conn.close()
    
    return None
    


# In[5]:

def get_pmids_and_article_data(query, db_name):
    '''Does same as get_pmids and get_article_data functions but in one step, which I think is faster than doing both separately. Uses the core resulttype from ePMC search module.
    Does not exclude articles that are in chembl because that field is set later.
    Return the query_id assigned to the query in the queries table, can then be used in subsequent functions.
    kwargs:
            query -- string
            db_name -- name of SQLite db
    
    '''
    
    moment = datetime.datetime.now()
    current_date = "{}-{}-{}".format(moment.year, moment.month, moment.day)
    
    base = 'http://www.ebi.ac.uk/europepmc/webservices/rest/search/query=({} AND (src:MED OR src:PMC OR src:CTX))&resulttype=core&page={}'
    
    conn = lite.connect(db_name)    
    cursor = conn.cursor()
    new_inserted_count = 0
    exist_count = 0
    
    def do_query(my_query, page_nr = 1):
        
        response = requests.get(base.format(my_query, page_nr))
        
        assert response.status_code == 200
        
        tree = etree.fromstring(response.content)
        
        return tree

    def retrieve_data(result):
    
        #pmid = result.xpath('pmid/text()')[0]
        title = result.xpath('title/text()')[0]
        abstract = result.xpath('abstractText/text()')[0]
        year = result.xpath('journalInfo/yearOfPublication/text()')[0]
        journal_title = result.xpath('journalInfo/journal/title/text()')[0]
        journal_abbrev_title = result.xpath('journalInfo/journal/medlineAbbreviation/text()')[0]
        
        if result.xpath('inEPMC/text()')[0] == 'Y': 
            in_ePMC = 1 
        else: 
            in_ePMC = 0

        url_list = []
        avail_code_list = []
        doc_style_list = []
        other_links_list = []
        pdf_links_list = []

        ft_url_list = result.xpath('fullTextUrlList/fullTextUrl')
        
        for url_item in ft_url_list:
            avail_code = url_item.xpath('availabilityCode/text()')[0]
            avail_code_list.append(avail_code)
            doc_style = url_item.xpath('documentStyle/text()')[0]
            doc_style_list.append(doc_style)

            if avail_code != 'S' and doc_style != 'pdf':
                url = url_item.xpath('url/text()')[0]
                other_links_list.append(url)

        if 'pdf' in doc_style_list:
            for url_item in ft_url_list:
                doc_style = url_item.xpath('documentStyle/text()')[0]
                if doc_style == 'pdf':
                    url = url_item.xpath('url/text()')[0]
                    pdf_links_list.append(url)                        
        else:
            pass

        avail_codes = ', '.join(set(avail_code_list)) if avail_code_list else None
        other_links = ', '.join(other_links_list) if other_links_list else None
        pdf_links = ', '.join(pdf_links_list) if pdf_links_list else None
        
        #cursor.execute("insert or ignore into result_ids(query_id, pmid) values(?,?)", (current_query_id, pmid))
        
        cursor.execute("insert or ignore into article_data(pmid, year, title, abstract, journal_title, journal_abbrev_title, in_epmc, avail_codes, pdf_links, other_links) values (?,?,?,?,?,?,?,?,?,?)", (pmid, year, title, abstract, journal_title, journal_abbrev_title, in_ePMC, avail_codes, pdf_links, other_links))
        
        conn.commit()
        
        return None

    tree = do_query(query)
 
    result_hitcount = int(tree.xpath('/responseWrapper/hitCount/text()')[0]) 
    cursor.execute("insert into queries(query_id, query, hitcount, date_performed) values (NULL,?,?,?)", (query, result_hitcount, current_date))
    current_query_id = cursor.lastrowid

    result_abstracts = tree.xpath('/responseWrapper/resultList/result')
    
    cursor.execute('select distinct pmid from article_data')
    existing_pmid_list = [i[0] for i in cursor.fetchall()]
    #print(len(existing_pmid_list))
    #print(existing_pmid_list[:10])

    for item in result_abstracts:
        try:
            pmid = int(item.xpath('pmid/text()')[0])
            #print(pmid)
            cursor.execute("insert or ignore into result_ids(query_id, pmid) values(?,?)", (current_query_id, pmid))

            if (pmid in existing_pmid_list) == True:
                #exist_count += 1
                continue
            else:
                try:
                    retrieve_data(item)
                    #new_inserted_count += 1
                except IndexError:
                    object_id = item.xpath('id/text()')[0]
                    error_comment = '(get_pmids_and_article_data) - IndexError with XML, possibly field not present, e.g. no journal info e.g. when is book chapter'
                    cursor.execute("insert or ignore into error_records(query_id, object_id, error_comment) values (?,?,?)", (current_query_id, object_id, error_comment))
                    continue
        except IndexError:
            object_id = item.xpath('id/text()')[0]
            error_comment = '(get_pmids_and_article_data) - IndexError with XML, possibly no pmid for this item or no journal info e.g. when is book chapter'
            cursor.execute("insert or ignore into error_records(query_id, object_id, error_comment) values (?,?,?)", (current_query_id, object_id, error_comment))
            continue
        
    if (result_hitcount % 25 == 0) and (result_hitcount / 25 > 1): # hitCount is multiple of 25 (remainder==0) but more than one page
        total_pages = int(result_hitcount / 25)

    else:  
        total_pages = (int(result_hitcount // 25)) + 1

    for page in range(2,total_pages+1): #page one was already done so continue from page 2 and add results to same dictionary

        tree = do_query(query, page_nr = page)
        result_abstracts = tree.xpath('/responseWrapper/resultList/result')

        for item in result_abstracts:
            try:
                pmid = int(item.xpath('pmid/text()')[0])
                cursor.execute("insert or ignore into result_ids(query_id, pmid) values(?,?)", (current_query_id, pmid))

                if pmid in existing_pmid_list:
                    #exist_count += 1
                    continue
                else:
                    try:
                        retrieve_data(item)
                        #new_inserted_count += 1
                    except IndexError:
                        object_id = item.xpath('id/text()')[0]
                        error_comment = '(get_pmids_and_article_data) - IndexError with XML, possibly field not present, e.g. no journal info e.g. when is book chapter'
                        cursor.execute("insert or ignore into error_records(query_id, object_id, error_comment) values (?,?,?)", (current_query_id, object_id, error_comment))
                        continue
            except IndexError:
                object_id = item.xpath('id/text()')[0]
                error_comment = '(get_pmids_and_article_data) - IndexError with XML, possibly no pmid for this item or no journal info e.g. when is book chapter'
                cursor.execute("insert or ignore into error_records(query_id, object_id, error_comment) values (?,?,?)", (current_query_id, object_id, error_comment))
                continue
            
    conn.commit()
    conn.close()
    #print(new_inserted_count, exist_count)
    
    return current_query_id


# In[9]:

def get_scores(query_id, db_name):
    '''For given query_id rank corresponding titles and abstracts using ChEMBL HeCaToS webservice (ChEMBL-likeness score). 
    Save in scores table. Excludes any pmids for which there is already a score in the scores table.
    If no score could be retrieved, is saved in error table together with pmid and query_id
    kwargs: query_id -- query_id from queries table, db_name'''
    
    title_abstract_dict = {}
    URL = 'http://scitegic.windows.ebi.ac.uk:9955/rest/HeCaToS_ChEMBLLIKE/{}/{}'
    score_dict = {}
    #my_count = 0

    conn = lite.connect(db_name)
    cursor = conn.cursor()
    
    cursor.execute('''select pmid, title, abstract from article_data 
                    where title is not null 
                    and abstract is not null 
                    and pmid not in (select distinct pmid from scores) 
                    and pmid in (select pmid from result_ids where query_id = ?) ''', (query_id,))
    #this sql query excludes things that already have a score
    
    results = cursor.fetchall()
    
    for item in results:
        pmid = item[0]
        title = item[1]
        abstract = item[2]
        title_abstract_dict[pmid]= [title, abstract]
    
    for article_pmid in title_abstract_dict:
        
        try:
            title = title_abstract_dict[article_pmid][0]
            abstract = title_abstract_dict[article_pmid][1]
    
            pre_title = str(base64.b64encode(title.encode('utf-8')))
            post_title = re.search("^b'(.*)'$", pre_title).group(1)
    
            pre_abstract = str(base64.b64encode(abstract.encode('utf-8')))
            post_abstract = re.search("^b'(.*)'$", pre_abstract).group(1)
    
            response_0 = requests.get(URL.format(post_title, post_abstract))
            response = response_0.json()  
            
            cl_score = float(response['score'])
            cursor.execute("insert or ignore into scores(pmid, score) values(?,?)", (article_pmid, cl_score))
            #my_count += 1
            
        except ValueError:
            error_info = '(get scores) error, status code: '+str(response_0.status_code)
            cursor.execute("insert or ignore into error_records(query_id, pmid, error_comment) values (?,?,?)", (query_id, article_pmid, error_info))
        except AttributeError:
            error_info = '(get scores) AttributeError'
            cursor.execute("insert or ignore into error_records(query_id, pmid, error_comment) values (?,?,?)", (query_id, article_pmid, error_info))
    
    conn.commit()
    conn.close()
    #print(my_count)
    return None
    
    


# In[21]:

def set_chembl_values(query_id, db_name):
    '''This function updates the article_data table and sets the in_chembl field by comparing pmid with chembl_pmids table. Needs to be done for each query_id for the sorting in the get_df to work well.
    kwargs:
            query_id
            db_name'''
    
    conn = lite.connect(db_name)
    cursor = conn.cursor()
    
    cursor.execute('update article_data set in_chembl = 1 where pmid in (select pmid from result_ids where query_id = ?) and pmid in (select distinct pmid from chembl_pmids)', (query_id,))
    cursor.execute('update article_data set in_chembl = 0 where pmid in (select pmid from result_ids where query_id = ?) and pmid not in (select distinct pmid from chembl_pmids)', (query_id,))
    
    conn.commit()
    conn.close()

    return None


# In[11]:

def get_availabilities(query_id, db_name):
    '''Want to check access to paper via campus subscriptions for those without 'F' or 'OA' in the availability_codes OR those without any availability codes, of a given query. 
    Check whether the SFX resolver has link to Full Text available, if so, save link in campus_links field of the article_links table, else, save link to 'request document' form in request_access field in article_links table.
    Excludes articles already in ChEMBL. Also skips pubmed ids which have article data already retrieved.
    kwargs: query_id -- query_id from queries table in queries_db
            db_name -- name of SQLite database'''
    
    conn = lite.connect(db_name)
    cursor = conn.cursor()
    avail_dict = {}
            
    def f(form):
    
        params = {x.attrib['name']: x.attrib['value'] for x in form.xpath('.//input[@type="hidden"]')}
    
        response = requests.get('http://wtgcsfx.hosted.exlibrisgroup.com/wtsc/cgi/core/sfxresolver.cgi', params=params)
    
        return response.status_code, response.url    
    
    
    sql =  '''select pmid from article_data
            where ((avail_codes not like '%F%' and avail_codes not like '%OA%') or avail_codes is null)
            and pmid in (select pmid from result_ids where query_id = ?)
            and pmid not in (select distinct pmid from article_links)
            and in_chembl != 1'''
    
    cursor.execute(sql, (query_id,))
   
    for pmid in [i[0] for i in cursor.fetchall()]:
        
        try:
        
            for attempt_number in range(3):

                response = requests.get('http://wtgcsfx.hosted.exlibrisgroup.com/wtsc?sid=Entrez:PubMed&id=pmid:{}'.format(pmid))

                if response.status_code == 200:
                    tree = lxml.html.fromstring(response.text)
                    ft_avail = tree.xpath('//div[@class="service"]/text()')
                    break

                else:
                    sleep(2)
                    continue


            if response.status_code == 200:

                if 'No Full text available' in str(ft_avail):

                    try:
                        #raise requests.TooManyRedirects('oh no')
                        responses = [f(x) for x in tree.xpath('//table[@id="service_type_header_getDocumentDelivery"]//form[contains(@name, "basic")]')]

                        urls = [response_tuple[1] for response_tuple in responses if response_tuple[0] == 200]

                        if not urls:
                            error_info = '(get_availability) Case 1: No URL with status_code = 200 available.'
                            cursor.execute("insert or ignore into error_records(query_id, pmid, error_comment) values (?,?,?)", (query_id, pmid, error_info))
                            continue

                        else:
                            urls_string = ', '.join(urls)
                            cursor.execute('insert into article_links(pmid, request_access) values(?,?)', (pmid, urls_string))
                            conn.commit()
                            continue

                    except requests.TooManyRedirects:
                        print(pmid, ' -- TooManyRedirects')
                        error_info = '(get_availability) TooManyRedirects'
                        cursor.execute("insert or ignore into error_records(query_id, pmid, error_comment) values (?,?,?)", (query_id, pmid, error_info))
                        conn.commit()
                        continue

                if not ft_avail:
                    error_info = '(get_availability) Case 2: ft_avail is empty. Not the usual SFX page, possibly SFX Multiple Object Menu, suggest to go to the page manually'
                    cursor.execute("insert or ignore into error_records(query_id, pmid, error_comment) values(?,?,?)", (query_id, pmid, error_info))
                    continue

                if ('No Full text available' not in str(ft_avail) and 'Request document via' in str(ft_avail)):

                    try:
                        #raise requests.TooManyRedirects('oh no')
                        responses = [f(x) for x in tree.xpath('//table[@id="service_type_header_getFullTxt"]//form[contains(@name, "basic")]')]

                        urls = [response_tuple[1] for response_tuple in responses if response_tuple[0] == 200]

                        if not urls:

                            try:
                                responses = [f(x) for x in tree.xpath('//table[@id="service_type_header_getDocumentDelivery"]//form[contains(@name, "basic")]')]
                                urls = [response_tuple[1] for response_tuple in responses if response_tuple[0] == 200]

                                if not urls:
                                    error_info = '(get_availability) Case 1b: No URL with status_code = 200 available.'
                                    cursor.execute("insert or ignore into error_records(query_id, pmid, error_comment) values (?,?,?)", (query_id, pmid, error_info))
                                    continue

                                else:
                                    urls_string = ', '.join(urls)
                                    cursor.execute('insert into article_links(request_access, pmid) values(?,?)', (urls_string, pmid))
                                    conn.commit()
                                    continue

                            except requests.TooManyRedirects:
                                #print(pmid, ' -- TooManyRedirects')
                                error_info = '(get_availability) TooManyRedirects'
                                cursor.execute("insert or ignore into error_records(query_id, pmid, error_comment) values (?,?,?)", (query_id, pmid, error_info))
                                conn.commit()
                                continue

                            #error_info = '(get_availability) Case 3: No campus ft-URL with status_code = 200 available.'
                            #cursor.execute("insert or ignore into other_error(query_id, pmid, error_comment) values (?,?,?)", (query_id, pmid, error_info))
                            #continue

                        else:
                            urls_string = ', '.join(urls)
                            cursor.execute('insert into article_links(campus_links, pmid) values (?,?)', (urls_string, pmid))
                            conn.commit()
                            continue

                    except requests.TooManyRedirects:
                        error_info = '(get_availability) TooManyRedirects'
                        cursor.execute("insert or ignore into error_records(query_id, pmid, error_comment) values (?,?,?)", (query_id, pmid, error_info))
                        conn.commit()
                        continues


                else:
                    error_info = "(get_availability) Case 4: Could not resolve access for article strings that were tested were not present in response."
                    cursor.execute('insert or ignore into error_records(query_id, pmid, error_comment) values (?,?,?)', (query_id, pmid, error_info))
                    continue

            else:
                error_info = "(get_availability) Case 5: could not get status code 200 from sfx page for this pmid after 3 attempts."
                cursor.execute('insert or ignore into error_records(query_id, pmid, error_comment) values (?,?,?)', (query_id, pmid, error_info))
                continue
        
        except (TypeError, requests.ConnectionError):
            print('oops, not going well there', pmid)
            continue
        except:
            print("Unexpected error:", sys.exc_info()[0], pmid)
            continue
            
    conn.commit()
    conn.close()
            
    return None
    
    
    
    
    


# In[20]:

def get_df(query_id_list, db_name, sql_condition = None):
    '''Makes a dataframe for a given query/queries. Can include multiple query_ids from queries table. If only one is needed put that one item in a list. 
    Selects following information on results for a given query from the queries_db: pmid, year, title, abstract, in_chembl, score, availability_codes, pdf_links, campus_links, request_access.
    Then sort the dataframe. First sorts on in_chembl, availability (subscription without access last), then on scores. Campus_links are only displayed if full text available from subscription.
    At the end of the dataframe the abstracts that had an error are displayed with an error_comment. The error_comment lists the function when the error occurred.
    Returns a tuple with first a non-HTML data frame and second an HTML version of the dataframe. Save the non-HTML version for use in the colour_terms function.
    It is possible to add further condition(s) to the sql statement by using the sql_condition argument. One "and" will be inserted by the function, then can add e.g. 's.score > 10', 
    which will then be appended to the sql statement.
    kwargs: query_id_list - list of query_ids to be included.
            db_name -- name of SQLite database
            sql_condition -- a further condition to be appended to the sql statement. One 'and' will be included by the function, so write the condition straight away. Default = None'''
        
    conn = lite.connect(db_name)
    cursor = conn.cursor()
    pd.set_option('max_colwidth',100000)
    
    sql = '''select distinct r.pmid, a.year, a.title, a.abstract, a.journal_title, a.in_chembl, s.score, a.avail_codes, a.pdf_links, a.other_links, al.campus_links, al.request_access, er.error_comment
        from result_ids r
        left join article_data a on r.pmid = a.pmid
        left join scores s on a.pmid = s.pmid
        left join article_links al on a.pmid = al.pmid
        left join error_records er on (r.query_id = er.query_id and r.pmid = er.pmid)
        where r.query_id in ({})
        '''
    
    if sql_condition != None:
        sql = sql+' and '+sql_condition
    
    
    def make_into_links(x, link_type):
        if x:
            link_list = x.split(', ')
            new_link_list = []
            if link_list:
                for index,i in enumerate(link_list):
                    item = '<a target="_blank" href="{}">{}_link_{}</a>'.format(i, link_type, index)
                    new_link_list.append(item)
            return ', '.join(new_link_list)
    
    sql_result_list = []
    
    query_id = str(query_id_list).strip('[]')
        
    for row in cursor.execute(sql.format(query_id)):
        sql_result_list.append(row)
    
                              
    try:
        df = pd.DataFrame(sql_result_list, columns = ['pmid', 'year', 'title', 'abstract', 'journal', 'inChEMBL', 'score', 'avail', 'pdf_links', 'other_links', 'campus_links', 'request_access', 'error_comment'])

        df['ePMC_link'] = df['pmid'].apply(lambda x: '<a target="_blank" href="http://europepmc.org/abstract/MED/{}">ePMC link</a>'.format(x))
        
        df['pdf_links'] = df['pdf_links'].apply(make_into_links, link_type = 'pdf')
        df['pdf_links_boolean'] = df['pdf_links'].apply(lambda x: 1 if x else 0)
        df['pdf_links'] = df['pdf_links'].apply(lambda x: x if x else '')
        
        df['other_links'] = df['other_links'].apply(make_into_links, link_type = 'other')
        df['other_links_boolean'] = df['other_links'].apply(lambda x: 1 if x else 0)
        df['other_links'] = df['other_links'].apply(lambda x: x if x else '')
        
        df['campus_links'] = df['campus_links'].apply(make_into_links, link_type = 'campus')
        df['campus_links_boolean'] = df['campus_links'].apply(lambda x: 1 if x else 0)
        df['campus_links'] = df['campus_links'].apply(lambda x: x if x else '')
        
        df['request_access'] = df['request_access'].apply(make_into_links, link_type = 'request_access')
        df['request_access_boolean'] = df['request_access'].apply(lambda x: 1 if x else 0)
        df['request_access'] = df['request_access'].apply(lambda x: x if x else '')
        
        df['avail'] = df['avail'].apply(lambda x: x if x else '')
        
        df['error_comment'] = df['error_comment'].apply(lambda x: x if x else '')
        
        df['pmid_link'] = df['pmid'].apply(lambda x: '<a target="_blank" href="http://www.ncbi.nlm.nih.gov/pubmed/?term={}">{}</a>'.format(x,x))
                                              
        #df = df.convert_objects(convert_numeric=True)

        df.sort(columns = ['inChEMBL', 'error_comment', 'request_access_boolean', 'score', 'pdf_links_boolean', 'other_links_boolean', 'campus_links_boolean'], ascending = [True, True, True, False, False, False, False], inplace = True)
        df.index = range(1,len(df) + 1)
        
        df = df.loc[:,['pmid', 'pmid_link', 'year', 'title', 'abstract', 'journal', 'inChEMBL', 'score', 'avail', 'pdf_links', 'other_links', 'campus_links', 'request_access', 'error_comment']]

        conn.close()
    
        return df, HTML(df.to_html(escape=False))

    except ValueError as e:
        print(e)
        print('probably no results for this query, empty db table or problem with df e.g. sorting')

    


# In[2]:

def colour_terms(df, markup_list):
    '''Will colour the terms in the df according to the markup_list and return an HTML object with the colours. Supply a df(non-HTML) and list of dictionaries. Each dictionary should have keys 'name'(of dict), 'terms'(value is list of terms), and 'colour'.
    kwargs: df -- a dataframe with a column 'abstracts' and 'title'
            markup_list -- a list of dictionaries with markup specification (list of terms and colour assigned)
    '''
    df = df.copy()
    
    css = '\n'.join(".{name} {{ color: {color}; }}".format(**x) for x in markup_list)

    def add_colour(html):
        
        for x in markup_list:

            pattern = re.compile(r'\b(' + '|'.join(x['terms']) + r')\b', flags=re.I)

            html = re.sub(pattern, r'<span class="{}">\1</span>'.format(x['name']), html)
            
        return html

    df['abstract'] = df['abstract'].apply(add_colour)
    df['title'] = df['title'].apply(add_colour)
    
    return HTML('<style>' + css + '</style>' + df.to_html(escape=False))


# In[22]:

def separate_column_df(query_id_list, db_name, sql_condition = None):
    '''Same as get_df function only now pdf_links and other_links have been split up in separate column for each link.
    Get_df docstring: Selects following information on results for a given query from the queries_db: pmid, year, title, abstract, in_chembl, score, availability_codes, pdf_links, campus_links, request_access.
    Then sort the dataframe. First sorts on in_chembl, availability (subscription without access last), then on scores. Campus_links are only displayed if full text available from subscription.
    At the end of the dataframe the abstracts that had an error are displayed with an error_comment. The error_comment lists the function when the error occurred.
    It is possible to add further condition(s) to the sql statement by using the sql_condition argument. One "and" will be inserted by the function, then can add e.g. 's.score > 10', 
    which will then be appended to the sql statement.
    Returns a tuple with first a non_HTML data frame and second an HTML version of the dataframe.
    kwargs: query_id_list -- query_id from queries table in queries_db
            db_name -- Name of SQLite database
            sql_condition -- a further condition to be appended to the sql statement. One 'and' will be included by the function, so write the condition straight away. Default = None'''
    
    conn = lite.connect(db_name)
    cursor = conn.cursor()
    pd.set_option('max_colwidth',100000)
    
    sql = '''select distinct r.pmid, a.year, a.title, a.abstract, a.journal_title, a.in_chembl, s.score, a.avail_codes, a.pdf_links, a.other_links, al.campus_links, al.request_access, er.error_comment
        from result_ids r
        left join article_data a on r.pmid = a.pmid
        left join scores s on a.pmid = s.pmid
        left join article_links al on a.pmid = al.pmid
        left join error_records er on (r.query_id = er.query_id and r.pmid = er.pmid)
        where r.query_id in ({})
        '''
    
    if sql_condition != None:
        sql = sql+' and '+sql_condition
    
    def make_into_links(x, link_type):
        if x:
            link_list = x.split(', ')
            new_link_list = []
            if link_list:
                for index,i in enumerate(link_list):
                    item = '<a target="_blank" href="{}">{}_link_{}</a>'.format(i, link_type, index)
                    new_link_list.append(item)
            return ', '.join(new_link_list)
    
    sql_result_list = []
    
    query_id = str(query_id_list).strip('[]')
    
    for row in cursor.execute(sql.format(query_id)):
        sql_result_list.append(row)
    
    try:
        df = pd.DataFrame(sql_result_list, columns = ['pmid', 'year', 'title', 'abstract', 'journal', 'inChEMBL', 'score', 'avail', 'pdf_links', 'other_links', 'campus_links', 'request_access', 'error_comment'])

        df['ePMC_link'] = df['pmid'].apply(lambda x: '<a target="_blank" href="http://europepmc.org/abstract/MED/{}">ePMC link</a>'.format(x))
        
        df['pdf_links'] = df['pdf_links'].apply(make_into_links, link_type = 'pdf')
        df['pdf_links_boolean'] = df['pdf_links'].apply(lambda x: 1 if x else 0)
        df['pdf_links'] = df['pdf_links'].apply(lambda x: x if x else '')
        
        df['other_links'] = df['other_links'].apply(make_into_links, link_type = 'other')
        df['other_links_boolean'] = df['other_links'].apply(lambda x: 1 if x else 0)
        df['other_links'] = df['other_links'].apply(lambda x: x if x else '')
        
        df['campus_links'] = df['campus_links'].apply(make_into_links, link_type = 'campus')
        df['campus_links_boolean'] = df['campus_links'].apply(lambda x: 1 if x else 0)
        df['campus_links'] = df['campus_links'].apply(lambda x: x if x else '')
        
        df['request_access'] = df['request_access'].apply(make_into_links, link_type = 'request_access')
        df['request_access_boolean'] = df['request_access'].apply(lambda x: 1 if x else 0)
        df['request_access'] = df['request_access'].apply(lambda x: x if x else '')
        
        df['avail'] = df['avail'].apply(lambda x: x if x else '')
        
        df['error_comment'] = df['error_comment'].apply(lambda x: x if x else '')
        
        df['pmid'] = df['pmid'].apply(lambda x: '<a target="_blank" href="http://www.ncbi.nlm.nih.gov/pubmed/?term={}">{}</a>'.format(x,x))
                                              
        #df = df.convert_objects(convert_numeric=True)

        df.sort(columns = ['inChEMBL', 'error_comment', 'request_access_boolean', 'score', 'pdf_links_boolean', 'other_links_boolean', 'campus_links_boolean'], ascending = [True, True, True, False, False, False, False], inplace = True)
        df.index = range(1,len(df) + 1)
        
        # get highest number of pdf links in the table
        pdf_links_lengths_list = []
        df['pdf_links'].apply(lambda x: pdf_links_lengths_list.append(len(x.split(', '))))

        # insert each pdf link into separate column
        for number in range(0,max(pdf_links_lengths_list)):
            column_name = 'pdf-'+str(number)
            df = pd.concat([df, df['pdf_links'].apply(lambda x: pd.Series({column_name:x.split(', ')[number] if (len(x.split(', ')) > number) else ''}))], axis = 1)
    
        # get names of columns that were added to the table
        pdf_link_columns = []
        for i in df.columns:
            if 'pdf-' in i:
                pdf_link_columns.append(i)
            
        # repeat process done on pdf_links on other_links
        
        
        # get highest number of pdf links in the table
        other_links_lengths_list = []
        df['other_links'].apply(lambda x: other_links_lengths_list.append(len(x.split(', '))))

        # insert each other link into separate column
        for number in range(0,max(other_links_lengths_list)):
            column_name = 'other-'+str(number)
            df = pd.concat([df, df['other_links'].apply(lambda x: pd.Series({column_name:x.split(', ')[number] if (len(x.split(', ')) > number) else ''}))], axis = 1)
    
        # get names of columns that were added to the table
        other_link_columns = []
        for i in df.columns:
            if 'other-' in i:
                other_link_columns.append(i)
        
        
        #subset dataframe to get relevant columns only
        
        column_list = ['pmid', 'year', 'title', 'abstract', 'journal', 'inChEMBL', 'score', 'avail', 'pdf_links_boolean']
        
        for item in pdf_link_columns:
            column_list.append(item)
        
        for item in other_link_columns:
            column_list.append(item)
            
        for item in "campus_links, request_access, error_comment".split(', '):
            column_list.append(item)
        
        print(column_list)
        df = df.loc[:,column_list]

        conn.close()
    
        return df, HTML(df.to_html(escape=False))

    except ValueError as e:
        print(e)
        print('probably no results for this query, empty db table or problem with df e.g. sorting')


# In[2]:

def plot_scores(query_id_list, db_name, figure_title):
    '''Creates a histogram of the chembl-likeness-scores for a given query/queries. Will plot inline.
    kwargs: query_id_list
            db_name
            figure_title'''
    
    get_ipython().magic('matplotlib inline')
    
    conn = lite.connect(db_name)
    cursor = conn.cursor()
    
    query_id = str(query_id_list).strip('[]')
    cursor.execute('select score from scores where pmid in (select distinct pmid from result_ids where query_id in ({}))'.format(query_id))
    
    score_list = [i[0] for i in cursor.fetchall()]
    score_array = np.array(score_list)
    
    plt.hist(score_array)
    plt.xlabel('ChEMBL-likeness score')
    plt.ylabel('Number of abstracts')
    plt.title('Query: {}\n'.format(figure_title)+'Total number of distinct abstracts ranked: '+str(len(score_list)))
    
    conn.close()
    


# In[ ]:



