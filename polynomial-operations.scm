#lang racket

(provide mono-divides?
		 order
		 leading-term
		 first-nonzero-coefficient
         pad-left
		 pad-right
		 zeros
		 make-monomial
		 leading-term
	     poly-add
		 poly-sub
		 mono-poly-mult
		 mono-divide
		 poly-mult
		 poly-exp
		 poly-quotient-remainder
		 poly-quotient
		 poly-remainder
		 GCD-poly)

;basic methods--------------------------
(define (accumulate init op sequence)
	(if (null? sequence)
		init
		(op (car sequence) (accumulate init op (cdr sequence)))))

(define (fold-left init op sequence)
	(define (iter result seq)
		(if (null? seq)
			result
			(iter (op result (car seq)) (cdr seq))))
	(iter init sequence))

(define (polynomial? f)
	(accumulate #t (lambda (x y) (and x y)) (map number? f)))

(define (mono-divides? a b)
	(and (not (= (order a) -1))
	     (<= (order a) (order b))))

(define (order f)
	(cond ((zeros? f)
		   -1)
		  ((not (= (car f) 0))
		   (- (length f) 1))
		  (else 
			(order (cdr f)))))
		
(define (coefficient m)
	;used for monomials
	(first-nonzero-coefficient m))

(define (first-nonzero-coefficient f)
	(if (not (= (car f) 0))
		(car f)
		(first-nonzero-coefficient (cdr f))))

(define (remove-front-zeros f)
	(if (or (= (length f) 1)
			(not (= (car f) 0)))
		f
		(remove-front-zeros (cdr f))))

(define (rest f)
	(if (= (length f) 0)
		f
		(cdr (remove-front-zeros f))))

(define (pad-left f n)
	(if (= n 0)
		f
		(pad-left (cons 0 f) (- n 1))))

(define (pad-right f n)
	(if (= n 0)
		f
		(pad-right (append f (list 0)) (- n 1))))

(define (zeros n)
	(if (= n 0)
		'()
		(cons 0 (zeros (- n 1)))))

(define (zeros? f)
	(if (null? f)
		#t
		(and (= (car f) 0) (zeros? (cdr f)))))

(define (make-monomial term o)
	(pad-right (list term) o))

(define (leading-term f)
	(make-monomial (first-nonzero-coefficient f) (order f)))

(define (same-length? f g)
	(= (length f) (length g)))

(define (make-same-length f g)
	;where the order of g is >= the order of f
	(pad-left f (abs (- (length f) (length g)))))

(define (fix-length f g procedure)
	(if (> (length f) (length g))
	    (procedure f (make-same-length g f))
	    (procedure (make-same-length f g) g)))
		
;----------------------------------------

;arithmetic methods----------------------

(define (poly-add f . g)
	(define (poly-add* a b)
		(if (not (same-length? a b))
			(fix-length a b poly-add)
			(map + a b)))
	(accumulate f poly-add* g))

(define (poly-sub f . g)
	(define (poly-sub* a b)
		(if (not (same-length? a b))
			(fix-length a b poly-sub)
			(map - a b)))
	(fold-left f poly-sub* g))

(define (mono-poly-mult m f)
	(pad-right (map (lambda (x) (* (coefficient (leading-term m)) x)) f)
			   (if (> (order m) -1)
					(order m)
					0)))

(define (poly-mult f . g)
	(define (poly-mult* a b)
		(if (or (null? a)
			    (= (order a) -1))
			(list 0)
		    (poly-add (mono-poly-mult (leading-term a) b)
					          (poly-mult* (rest a) b))))
	(accumulate f poly-mult* g))

(define (poly-exp f n)
	(cond ((= n 0)
			(list 1))
		  ((even? n)
			(poly-exp (poly-mult f f) (/ n 2)))
		  (else (poly-mult f (poly-exp f (- n 1))))))


(define (mono-divide m n)
	(make-monomial (/ (first-nonzero-coefficient m) (first-nonzero-coefficient n)) 
				   (- (order m) (order n))))

(define (poly-quotient-remainder f g)
	(define (iter q r)
		(if (or (= (order r) -1)
			 	(not (mono-divides? (leading-term g) (leading-term r))))
			(list q r)
			(iter (poly-add q (mono-divide (leading-term r)
										   (leading-term g)))
				  (poly-sub r (poly-mult (mono-divide (leading-term r) 
													  (leading-term g))
										  g)))))
	(iter (list 0) f))

(define (poly-remainder f g)
	(cadr (poly-quotient-remainder f g)))

(define (poly-quotient f g)
	(car (poly-quotient-remainder f g)))

(define (GCD-poly f . g)
	(define (GCD-poly* a b)
		(if (= (order b) -1)
			a
			(GCD-poly* b (poly-remainder a b))))
	(accumulate f GCD-poly* g))
;------------------------------------------
