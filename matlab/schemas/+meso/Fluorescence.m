%{
# fluorescence traces before spike extraction or filtering
-> meso.Segmentation
%}


classdef Fluorescence < dj.Computed

	methods(Access=protected)

		function makeTuples(self, key)
		%!!! compute missing fields for key here
% 			 self.insert(key)
		end
	end

end