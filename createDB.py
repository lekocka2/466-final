#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 20:54:36 2021

@author: logankocka
"""

import pymysql

db = pymysql.connect(host='localhost', 
                     user='kockale', 
                     password='bio466', 
                     database='kockale')
cur = db.cursor()

sql = "DROP TABLE IF EXISTS summary"
cur.execute(sql)
db.commit()

sql = "DROP TABLE IF EXISTS overlaps"
cur.execute(sql)
db.commit()

#### create tables

createSummary = ("create table summary (num_CG_pairs_sense int NOT NULL, num_CG_pairs_antisense int, total_num_islands int, avg_island_length")
cur.execute(createSummary)
db.commit()

createOverlaps = ("create table overlaps (geneID VARCHAR(17) NOT NULL, gene_start int, gene_end int, island_start int, island_end int")
cur.execute(createOverlaps)
db.commit()


