#!/bin/bash

gofem coarse && GenVtu coarse
gofem fineNonlin && GenVtu fineNonlin
