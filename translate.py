from copy import deepcopy
import params

def frame_translate(seq,frame):
	out_seq = deepcopy(seq)
	out_seq.id = 'Frame: '+ str(frame)
	out_seq.description = ""

	if frame>0:
		out_seq.seq = out_seq.seq[frame-1:].translate()
	else:
		out_seq.seq = out_seq.seq.reverse_complement()[abs(frame)-1:].translate()
	return out_seq

def translate(seq):
	seqs = [frame_translate(seq,frame) for frame in params.frames]
	return seqs
