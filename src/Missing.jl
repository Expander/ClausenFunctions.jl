for f in (:cl1, :cl2, :cl3, :cl4, :cl5, :cl6)
    @eval $(f)(::Missing) = missing
end

for f in (:cl, :sl)
    @eval $(f)(::Integer, ::Missing) = missing
end
