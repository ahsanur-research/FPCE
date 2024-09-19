all: fpce fpce_wo_eb grh2odes

#modified version of FPCE that turns off the APPLY_EDGE_BOUND feature
fpce_wo_eb: fpce-wo-edge-bound.c problem.c itemset.c queue.c sgraph.c stdlib2.c aheap.c alist.c undo.c base.c vec.c
	gcc -Wall -O3 -Os fpce-wo-edge-bound.c -lm -o fpce_wo_eb


grh2odes: grh2odes.c problem.c itemset.c queue.c sgraph.c stdlib2.c aheap.c alist.c undo.c base.c vec.c
	gcc -Wall -O3 -Os grh2odes.c -lm -o grh2odes

fpce: fpce.c problem.c itemset.c queue.c sgraph.c stdlib2.c aheap.c alist.c undo.c base.c vec.c
	gcc -Wall -O3 -Os fpce.c -lm -o fpce

clean:
	rm fpce fpce_wo_eb grh2odes fpce_w_ob


