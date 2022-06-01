function [X, inds] = roundtowardvec(X,roundvec)
    sz=size(X);
    %Calculate differences
    X=X(:);
    roundvec=roundvec(:)';
    diffs=bsxfun(@minus,X,roundvec);
    [~,inds]=min(abs(diffs),[],2);
    X=roundvec(inds);
    X=reshape(X(:),sz);
    inds=reshape(inds(:),sz);
end