#!/bin/bash
cd _drafts
date=`date +"%Y-%m-%d"`

mv $1 ../_posts/$date-$1
