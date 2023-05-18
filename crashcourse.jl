### A Pluto.jl notebook ###
# v0.19.20

using Markdown
using InteractiveUtils

# ╔═╡ a22dba1a-9881-11ed-3f52-559778073a9c
using  PlutoUI,Images

# ╔═╡ e029e588-7e17-4888-bb14-a7aa94a95f82
using Plots

# ╔═╡ 066a6395-81cb-4323-85da-31deccc9f7e0
using LaTeXStrings

# ╔═╡ df6f74e1-42eb-4e28-9b7f-1d8bb8bb5b0f
html"<button onclick='present()'>Presentación</button>"

# ╔═╡ 331513e1-d6c0-4f4f-819b-4e32319bce80
md"""
# Programación en Julia: Herramienta para la enseñanza
${\Huge\text{\bf Julia Crash Course}}$

**Benjamín Pérez**

**Héctor Medel**
"""

# ╔═╡ b11c332f-2375-45f7-a4ab-72ed58a5c3f8
md"""
## Empezando con VSCode
- Podemos abrir VSCode haciendo click sobre el ícono de la aplicación
"""

# ╔═╡ af48b15a-f04c-42ab-a1fe-8ff965e0a6d4
load("/Users/hector/Desktop/Captura de pantalla 2023-01-22 a la(s) 21.41.19.png")

# ╔═╡ e5a72dd2-2d29-42dc-b125-93e7072c2697
md"""
- Abrimos un nuevo archivo de Julia
"""

# ╔═╡ 802e3512-814d-48b7-b5b6-54ecbdcf63bb
load("/Users/hector/Desktop/Captura de pantalla 2023-01-22 a la(s) 21.43.32.png")

# ╔═╡ 95c959ab-9900-4290-98ac-57af9bec14a9
md" - Seleccione el lenguaje: Julia en este caso"

# ╔═╡ 754d000d-1f35-4a29-ad1e-0ab6316e33a1
load("/Users/hector/Desktop/Captura de pantalla 2023-01-22 a la(s) 21.49.54.png")

# ╔═╡ 0966daf4-8e14-4569-b68a-54f24dbf04e5
md"""
- Iniciar REPL 


```julia
Shift + Enter
```

"""

# ╔═╡ 706b8cc1-5b96-4bda-ae4f-fd32d6842187
load("/Users/hector/Desktop/Captura de pantalla 2023-01-22 a la(s) 22.08.34.png")

# ╔═╡ b612721a-1b06-49a8-86b7-8a3f09480597
md"""
### Algunos atajos

#### Sobre el editor
- Ejecutar código en REPL y moverse: ```Shift+Enter```
- Ejecutar código en REPL: ```Ctrl+Enter```
- Ejecutar código de una celda en REPL: ```Alt+Enter``` 
- Ejecutar código de una celda en REPL y moverse: ```Alt+Shift+Enter``` 
- Limpiar resultados de la línea Clear Inline Results In Editor: Ctrl+I Ctrl+C
#### Sobre Terminal
- Interrumpir Ejecución: Ctrl+C
- Limpiar los resultados de la línea actual: ```Escape```
- Buscar en el historial: ```Ctrl+r```
- Moverse en el historial: ```arriba o abajo```
- Salir: ```Ctrl+d```
"""

# ╔═╡ 0edc622e-cdd4-4595-a43b-eed80ab301a8
md"""
### Correr código de un script
Podemos escribir código de ```Julia``` en el editor. Los comandos se autocompletarán.

###### Ejemplo: Hola mundo
```julia
print("Hola mundo!")
```

"""

# ╔═╡ 5f1e7909-583b-42d1-a9bb-9e63ae7ceb99
md"""
- Salvemos el archivo con algún nombre, por ejemplo: ```Hola1.jl``` 
- Corramos el código en REPL
```julia
include("Hola1.jl")
```
"""

# ╔═╡ ed620ff7-9458-4d8e-9d51-62fe0a97c623
md"""
### Símbolos matemáticos
Podemos usar símbolos de LaTeX
###### Ejemplo: Símbolos griegos
```julia
θ = π/3
```
```
\theta+Tab
\pi+Tab
```
"""

# ╔═╡ 6d2c74a1-3fd1-4310-b6a1-c2eac4243257
md"""
## Variables y tipos
Julia es un lenguaje de tipado opcional; se puede escribir un código sin mencionar el tipo de variable.
"""

# ╔═╡ 00fe7ab2-b468-4489-aade-9935d377d236
md"""
###### Ejemplo:

```julia
a = 1 + 1;
typeof(a) # Int64

a = "hola"
typeof(a) # String
```

"""

# ╔═╡ 83fa2bd9-b168-46c2-87d2-2d4fa2844e47
md"""
### Jerarquía de tipos
"""

# ╔═╡ 0adba4a5-c418-43b5-aa08-6b29aaf025c1
load("/Users/hector/Desktop/Julia-number-type-hierarchy.svg.png")

# ╔═╡ 06f24748-40fb-480f-a202-4de665de2ec8
md"""
### Conversion de enteros a flotantes
```julia
Float64(2) # 2.0 (double-precision)
Float32(2) # 2.0f0 (single-precision)
Float16(2) # Float16(2.0) (half-precision)
```

"""

# ╔═╡ 65bc4440-574c-409b-95ce-970584670d42
md"""
### Conversión de flotante a entero
```julia
Int64(2.0)        # 2
Int64(2.4)        # Error!
floot(Int64,2.4)  # 2
ceil(Int64,2.4)   # 3
round(Int64,2.4)  # 2
```

"""

# ╔═╡ 02d484f5-0090-42a8-8c89-7f76b19d4ac5
md"""
### División
Existen tres operaciones asociadas
- División normal ```/```
- División truncada a un entero ```div(x,y)```
- El resto de una división ```rem(x,y)```

###### Ejemplos
```julia
a = 1/2  # 0.5 (Float 64)

div(10,3) # 3
÷(10,3)   # 3

rem(10,3)   # 1
10%3        #1
```

"""

# ╔═╡ d8ee737b-3e56-4f71-9143-cbef1c4ec618
md"""
### Variables Booleanas
Las variables boolenas son ```true``` y ```false```.

#### Operadores lógicos
Los operadores lógicos son



| Operador | Sintáxis | Descripción |
|:---------|:---------|:------------|
| ```not``` |```!x```|regresa el opuesto de x|
| ```and``` | ```x&&y```| regresa ```true```si los dos operadores son ciertos |
| ```or```  | ```x \|\| y```| regresa ```true```si al menos uno de los operandos es cierto |




"""

# ╔═╡ 9741bd55-5d21-4a1e-9f0f-2819de092ae2
md"""
### Booleano a número
Como Bool es un subtipo de ```Entero```, entonces, ```true``` regresa ```1``` y ```false``` regresa ```0```

###### Ejemplo:
```julia
true + true # 2
```

"""

# ╔═╡ b3d8fcab-9413-43cd-a032-80f837dd6896
md"""
### Strings
Julia tiene un buen manejo de strings

###### Ejemplo: Strings y chars

```julia
"Hola"    # string
'H'       # char
'Hola'    # Error
```
"""

# ╔═╡ 82c606ec-06ea-4ab8-89e3-a84e314bc5bd
md""" 
#### Concatenación 

###### Ejmplo:
```julia
"Hola" * "mundo"            # "Holamundo"
string("Hola"," ","mundo")  # "Hola mundo"
```

"""

# ╔═╡ e929b471-563c-49db-8c6d-01b5924e4e03
md"""
#### Convertir un número a string
Hay varias formas de convertir un número a string
###### Ejemplo:
```julia
string(1/7)            #Regresa "0.14285714285714285"
"$(1/7)"               #Regresa "0.14285714285714285"


using Printf
@sprintf("%6.4f",1/7)  #Regresa "0.1429"
```
"""

# ╔═╡ 9be54d29-b85e-4337-bf0a-bb2cc1e13e02
md""" 
#### Función ```string```
###### Ejemplo
```julia
string("podemos concanetar una string a un número ",5/6) # "podemos concanetar una string a un número 0.8333333333333334"
```
"""

# ╔═╡ 82dfe60b-2779-4ec0-96eb-023b9a1d84d1
md"""
##### Interpolación
```julia
"podemos concanetar una string a un número $(5/6)" # "podemos concanetar una string a un número 0.8333333333333334"
```

"""

# ╔═╡ 63573df9-2133-4e10-89b6-b33a7fe5d3de
md"""
## Bucles
"""

# ╔═╡ f11eaf6a-e3d0-4ed5-b64e-1b0a27ffe398
md"""
### Bucle simple ```for```
"""

# ╔═╡ e3e66c98-c462-4cb9-86f0-d462a005977b
md"""
###### Ejemplo
```julia
x = 1:10;

for i in x
 println(i)
end
```
"""

# ╔═╡ 7f01e0e9-9885-4d85-802a-457980507d87
let
	x = 1:10;
	for i in x
		println(i)
	end
end

# ╔═╡ f8bca9af-a155-496d-a1c7-2eb0b0822037
md"""
Ejemplo:
```julia
x = 0;
for k in 1:1000
 x = x + (1/k)^2
end
```

"""

# ╔═╡ ed3951f4-df95-468f-94b9-c694d13c1605
let
	x = 0;
for k in 1:1000
 x = x + (1/k)^2
end
	println(x)
end

# ╔═╡ ffd77c4b-99c6-4591-95eb-15ca69ed5361
md"""
### Bucles anidados
```julia
for i in 1:3
	for j in 1:3
		println("i = ", i, " j = ",j, "\n")
	end
end
```
"""

# ╔═╡ eaadf6cf-0cf5-43d9-943d-3975fc8a356d
for i in 1:3
	for j in 1:3
		println("i = ", i, " j = ",j, "\n")
	end
end

# ╔═╡ 182e369d-b31b-4f7b-b455-55cdfc0fb4b1
md"""
Otra forma equivalente es:
```julia
for i in 1:3, j  in 1:3
	println("i = ", i, " j = ",j, "\n")
end
```
"""

# ╔═╡ 4a094840-92f5-4447-aa9d-42690f40dcb2
md"""
#### Comandos ```break``` y  ```continue```
Sirven para romper un bucle o para hacer ciertas acciones sin considerar algunos valores

###### Ejemplo 1
```julia
x=0
for k in 1:100000
    term = (1/k)^2
    x = x + term
    if (abs(term) < 1e-10) 
		break 
	end
end
```

"""

# ╔═╡ ee262cee-ab01-49cb-bc70-2c73af12f1cd
let
x=0
for k in 1:100000
    term = (1/k)^2
    x = x + term
    if (abs(term) < 1e-10) break end
end
	println(x)
end

# ╔═╡ ab3980f9-a467-4212-8c18-dd1b2442f07c
md"""
###### Ejemplo
Lo anterior también lo podemos escribir con un bucle condicionado ```while```

```julia
x=0
iter = 0
while ( iter == 0 || abs(term) < 1e-10) && (iter < 100000)
    term = (1/k)^2
    x = x + term
    iter = iter + 1
end
```

"""

# ╔═╡ c944d642-b67a-41fa-a193-c156fdc02479
md"""
#### `continue`

Si queremos saltar una (o varias) repeticiones dentro de un loop podemos usar el comando `continue`.

###### Ejemplo
```julia
	for n in 1:10
		if 3<= n <= 6
			continue
		end
		println(n)
	end
```

"""

# ╔═╡ 87f969e4-e208-4dcf-b548-a1684345969e
let
		for n in 1:10
		if 3<= n <= 6
			continue
		end
		println(n)
	end
end

# ╔═╡ 02aa2202-5ec1-404e-a7bf-33b974c47588
md"""
## Funciones
"""

# ╔═╡ f39e45d1-33f3-4a5a-b152-5f39be2d6f89
md"""
#### Definiendo una función
Una función es un objeto que toma como *input* algunos argumentos, les hace operaciones (cualquiera permitida) y regresa valores.

Nota: Pueden ser de distintos tipos y estructuras. Los argumentos pueden estar expresados por tipo (opcional pero recomendado); los tipos pueden ser definidos por uno mismo.
"""

# ╔═╡ d8928380-e9d8-4004-a8d2-fae26cbe6a6e
md"""
La sintaxis general de una función es la siguiente:
```julia
function funcname(argumentos)
		#Something
		return values
end
```
"""

# ╔═╡ 05ee5b6d-b6f0-43eb-a794-ad5131d6b1f4
md"""
###### Ejemplo 1: Función zeta de Reimman
```julia
function sum_zeta(s,nterms)
    x = 0
    for n in 1:nterms
        x = x + (1/n)^s
    end
    return x
end
```
"""

# ╔═╡ c12fcbbd-46d3-413b-aff2-7631ae108cc4
md"""
```julia
sum_zeta(2,10000)
```
"""

# ╔═╡ d5b049b8-15ac-46a0-9145-bd3ea24c8aac
let
function sum_zeta(s,nterms)
    x = 0
    for n in 1:nterms
        x = x + (1/n)^s
    end
    return x
end	

	sum_zeta(2,10000)
end

# ╔═╡ 3c3b5816-928c-4c00-959d-6c33da5ec111
md"###### Ejemplo 2: Varios valores de salida"

# ╔═╡ 0360935c-09dc-4775-94bb-bc54dd631138
md"""
```julia
function multi(n,m)
	n*m, div(n,n), n%m
end
```
"""

# ╔═╡ 07f1cba5-eae0-4b5f-8c62-3557a71accc1
md"""
```julia
multi(8,2)
```
"""

# ╔═╡ 7d440a4f-2a47-453f-a47f-15a67964cc36
let
	function multi(n,m)
		n*m, div(n,n), n%m
	end

	multi(8,2)
end

# ╔═╡ a9a5b3f6-2019-4989-b04c-f55d2437f530
md"""
###### Ejemplo: Funciones que modifican su input (`!`)
Regularmente los valores no son copiados cuando pasan a una función. Una forma de decirle que lo puede cambiar es con el signo de exclamación `(!)`

```julia
function add_one!(x)
x .= x .+ 1
end
```

"""

# ╔═╡ 920bfd76-525f-4a5f-b4f5-2a330e620c8c
begin
function add_one!(x)
    x .= x .+ 1
end
	x = [1,2,3]
	add_one!(x)
end

# ╔═╡ 3fcee3a8-14ad-4dba-a6b6-b729ae0639e4
md"""
#### Funciones como expresiones matemáticas
```julia
f(x,y) = x^3 -x*y + 1/y;
f(4,5)
```

"""

# ╔═╡ 0d8b4d11-664e-4cd1-8d24-704870f65a1c
let
	f(x,y) = x^3 -x*y + 1/y;
	f(4,5)
end

# ╔═╡ 040e1ac8-b912-4e2b-a5de-8053133a0882
md"""
#### Funciones anónimas
- Se pueden definir funciones sin nombre


    function (x)
		x + 2
    end

- Funciones Lambda


    (x) -> x + 2 

- Funciones Lambda


    x -> x + 2
"""

# ╔═╡ 91e114d7-5347-462e-a06c-dd05fcb106c7
md"""
#### Funciones de funciones

Una función puede tomar una función como argumento
```julia
function ff(f::Function,x::Float)
	#Operaciones con f(x)
end
```
"""

# ╔═╡ a9745a65-c58b-49ab-b47d-a1df7466cad1
md"""
###### Ejemplo: Función de función
```julia
function derivada(f::Function,x::Float64,dx::Float64=0.001)
		df = (f(x+dx) - f(x-dx))/(2*dx)
		return df
end

f = x-> 2*x^2 + 30*x + 9; 

derivada(f,2.0,0.01)
```

"""

# ╔═╡ 2a79ca9b-dc5a-4198-990e-69967d0e3c82
let
function derivada(f::Function,x::Float64,dx::Float64=0.001)
		df = (f(x+dx) - f(x-dx))/(2*dx)
		return df
	end

	f = x-> 2*x^2 + 30*x + 9; 

	derivada(f,2.0,0.01)
end

# ╔═╡ 3b06956a-622f-424c-995c-679425c354c0
md"""
## Arreglos, vectores y matrices
"""

# ╔═╡ 71362aea-97eb-4fbf-b5f7-cb3644a70045
md"""
Las matrices tienen la siguiente forma:
```julia
A = [1 2 3; 1 2 4; 2 2 2]
```
"""

# ╔═╡ 662e3bb5-32ac-4f79-ba9a-51465b420a9e
let
	 A= [1 2 3; 1 2 4; 2 2 2]
end

# ╔═╡ fa4a472f-0972-46ba-b31c-51bd9ba23ff2
md"""
Es equivalente a
```julia
A = [1 2 3; 
     1 2 4; 
     2 2 2]
```

"""

# ╔═╡ d9cdaec4-91ac-4105-93f9-d6c393a3fa73
md"""
Otros arreglos se pueden escribir de la siguiente manera:
```julia
b1 = [4.0, 5, 6]                # 3-element Vector{Float64}
b2 = [4.0; 5; 6]                # 3-element Vector{Float64}
m1 = [4.0 5 6]                  # 1×3 Matrix{Float64}
```
"""

# ╔═╡ e9d1af7e-c8d4-4231-a600-373b56d7cefe
md"""
#### Arreglos `Any` 
Julia puede manejar arreglos con tipos no-numéricos

```julia
A = ["Hello", 1, 2, 3]
  4-element Vector{Any}:
  "Hello"
  1
  2
  3 
```
"""

# ╔═╡ 18479558-da15-4343-a51c-83aa20e92eed
md"""
#### Arreglos indefinidos
Se puede inicializar un arreglo vacío pero con las dimensiones bien definidas

```julia
n = 5
A1 = Array{Float64}(undef,n,n)          # 5×5 Matrix{Float64}
A2 = Matrix{Float64}(undef,n,n)         # 5×5 Matrix{Float64}

V1 = Array{Float64}(undef,n)            # 5-element Vector{Float64}
V2 = Vector{Float64}(undef,n)           # 5-element Vector{Float64}

A = Array{String}(undef,n)
A = Array{Any}(undef,n)
```

"""

# ╔═╡ 8eab62ad-4cc3-4b76-b02e-2f5b942be980
md"""
#### Arreglos vacíos
Podemos inicializar arreglos indefinidos y vacíos
```julia
v = Array{Float64}(undef,0)

v = Array{Float64}(undef,0)

v = []    #  Any[]
```

"""

# ╔═╡ f70276aa-3e9c-42f8-9021-46fccdf3802c
md"""
#### Funciones sobre arreglos y broadcasting (Operador `.`)
En Julia es muy sencillo aplicar funciones a arreglos.
"""

# ╔═╡ b8883753-a281-4c56-932e-b39987a67f97
md"###### Ejemplo: Operador ."

# ╔═╡ 1abf3f43-9f3d-407d-ab26-3aef2efe5d58
let
	f(x) = 3*x^3 / (1+x^2)
	x = [2π/n for n=1:30]
	y = f.(x)
end

# ╔═╡ 52cf0faf-4bff-4f55-88a1-e684810cdf03
md"También podemos escribirlo de la siguiente manera"

# ╔═╡ 34613b4d-34fa-4687-89f0-79d8a150486c
let
	y = @. 2*x^2 + 3*x^5 - 2*x^8
end

# ╔═╡ 24d42c25-3457-43c3-8ff6-d74a7719113a
md"""
### Indizado y rebanado
Es fácil obtener elementos de arreglos a través del indizado. Julia está muy optimizado para hacer ese tipo de labores.
"""

# ╔═╡ 678268a0-9dc2-4f5f-a3fb-6ee906401aef
md"###### Ejemplo: Primer y último elemento"

# ╔═╡ 1f2a56c3-f5de-453c-b965-9c066fd67aaa
let
	A = rand(6)
	b = A[begin]
	e = A[end]
	println(A)
	println()
	println("primer elemento es $b"," y ","último elemento es = $e")
end

# ╔═╡ b48aac04-5212-4f11-9884-b45564157174
md"###### Ejemplo: Rebanado"

# ╔═╡ 04950b12-3b15-4ac5-96f6-18963a3814a0
let
	A = rand(6,6)
	display(A)
	B = A[begin:2:end,begin:2:end]
	display(B)
	C = A[1:2:5,1:2:5]
	display(C)
end

# ╔═╡ 8ec133a1-76ac-4bca-924e-6880dfc6633c
md"###### Ejemplo: Indizado lógico"

# ╔═╡ 021e59a9-bd36-46e1-aa79-67e3a1953dfe
let
	A =  rand(6,6)
	display(A)
	A[A .< 0.5] .= 0
	display(A)
	
end

# ╔═╡ ff67f2aa-f208-48a9-8422-0d795145af62
md"###### Ejemplo: Iteraciones sobre arreglos"

# ╔═╡ bde8900d-6907-4b2c-b0f6-44428d83a684
let
	A = rand(6)
		for i ∈ eachindex(A)
    		println(string("i=$(i) A[i]=$(A[i])"))
		end
end

# ╔═╡ 11d479a7-71a2-4879-bc16-dc6053cb8c76
md"""
###### Otros comandos útiles
- `firstindex(A,dim)`
- `lastindex(A,dim)`
- `similar(Array{Float64}, axes(A))`
"""

# ╔═╡ fdfcd417-1d9c-4fd4-9df8-d50a59bddc1b
md"#### Operaciones con arreglos"

# ╔═╡ 5bd0a998-bcbd-4543-a899-480d0be845a7
md"###### Ejemplo: Multiplicación"

# ╔═╡ e51bd7c6-8649-4140-be0c-ae66902e834a
let
	A = rand(3,3)
	display(A)
	B = rand(3,3)
	display(B)
	C = A*B
	display(C)
end

# ╔═╡ 6faf62b9-003d-487e-98fe-23f5940c7562
let
	A = rand(4,4)
	display(A)
	v = rand(4)
	display(v)
	w = A*v
	display(w)
end

# ╔═╡ 7abe3e0f-6b6d-47a4-a86b-a7b666715ac2
md"###### Ejemplo: Multiplicación elemento por elemento"

# ╔═╡ ac884bda-f2a2-4666-bbd5-76339087dde3
let
	A = rand(3,3)
	display(A)
	B = rand(3,3)
	display(B)
	C = A.*B
	display(C)
end

# ╔═╡ 34ec3b4d-2d61-4666-ab76-9d71c274d582
md"###### Ejemplo: Producto punto"

# ╔═╡ 980a69af-e894-4d19-b5bf-9fac8f2900f5
let
	v = rand(100)
	w = rand(100)
	z = dot(v,w)
end

# ╔═╡ 27b65217-9841-4d6e-b60e-345b0c5e1245
md"###### Ejemplo: Operador backslash"

# ╔═╡ 25585d86-1909-47f7-92b8-297ef9b8c10b
let
	A = rand(3,3)
	b1 = [4.0, 5, 6]                # 3-element Vector{Float64}
	b2 = [4.0; 5; 6]                # 3-element Vector{Float64}
	m1 = [4.0 5 6]                  # 1×3 Matrix{Float64}

	x=A\b1                          # Solves A*x=b
	x=A\b2                          # Solves A*x=b  
	#x=A\m1                          # Error!!
end

# ╔═╡ b2b54e1f-89b9-4f0b-9d73-3d14dd6a7a2a
md"""
### Rearreglo y concatenación
"""

# ╔═╡ a65cb56f-0e37-47d3-939b-6b72aa23f230
md"""
#### Rearreglo en una dimensión
Podemos usar `push!` y  `pop!`
"""

# ╔═╡ c9b4affa-1e71-4b44-a92e-d914013e8747
md"###### Ejemplo"

# ╔═╡ f9119e93-da01-488b-817b-f8656e3b72b8
let
A = Float64[]       # Equivalent to A=Array{Float64}(undef,0)
push!(A, 4)         # Adds the number 4 at the end of the array
push!(A, 3)         # Adds the number 3 at the end of the array
v = pop!(A)         # Returns 3 and removes it from A	
end

# ╔═╡ 0cf416ca-ecc0-4193-9c03-18bdfb50d737
md"""
#### Concatenación
Podemos usar los comandos `hcat`, `vcat`, `cat`
"""

# ╔═╡ f4b3fec0-6786-415a-b709-9f311c50a6b2
md"""###### Ejemplo: `vcat`, `hcat` """

# ╔═╡ aca8a8ff-2708-4ccc-b6c7-027697bc5b7f
let
	A = [4 5 6]
	B = [6 7 8]

	M1 = vcat(A,B)
	display(M1)

	M2 = hcat(A,B)
	display(M2)
end

# ╔═╡ 94842940-9d38-4a36-bf26-382fa5ce15ef
md"""
###### Ejemplo: cat
La función `cat` es más general y permite realizar operaciones para cualquier dimensión
"""

# ╔═╡ 0e8ab455-8bdd-4688-b483-3e8cd3cb0b8f
let
	A = [4 5 6]
	B = [6 7 8]
	M1 = cat(A,B, dims=1)
	display(M1)

	M2 = cat(A,B, dims=2)
	display(M2)

	M3 = cat(A,B, dims=3)
	display(M3)
end

# ╔═╡ a22a5694-672b-489d-a272-b4aafa7c2d9d
md"#### Concatenación de matrices"

# ╔═╡ 89813b0d-1884-4888-b540-fa938d8ac0f4
md"Ejemplo"

# ╔═╡ 55d97f3b-873f-410d-b245-a4222072d6b6
let
	A = [4 5 6]
	B = [6 7 8]

	M1 = [A;B]
	display(M1)

	M2 = [A B]
	display(M2)
	
end

# ╔═╡ bc1e6fbc-4a1a-47cf-ae79-840fad31c1fb
md"""
## Estructura de datos
Algunas de las estructuras más usadas son: tuples, named tuples y dictionaries.
"""

# ╔═╡ 0628ca0c-1b8e-4f01-aa79-e7c99bc8f5f7
md"#### Inicializando tuplas"

# ╔═╡ 65449b04-9873-4c1d-aa8a-3a3133823424
md"###### Ejemplo:"

# ╔═╡ 26736f5a-64b8-4f39-a8cd-71218a922726
let
	t = (3.14,2.8) # Tuple{Float64,Float64}
	display(t)
end

# ╔═╡ fbf289d9-e3ca-43c1-96ad-e33f35addbf5
let
	t = 3.24, 2.8
	display(t)
end

# ╔═╡ 3a129775-7338-4c2b-b1de-92cf3ea5897e
md"#### Convertir tuplas a arreglos"

# ╔═╡ c1c33a9b-43be-4257-b44f-ada3bbfc3f50
md"###### Ejemplo"

# ╔═╡ 47be421e-24c0-453c-a8b6-eb9a118a26ec
let
	a = (1,2,3)
	t1 = collect(a)
	display(t1)

	t2 = [x for x in a]
	display(t2)

	t3 = [a...]
	display(t3)
end

# ╔═╡ 552daa2c-4ff9-4d03-b697-bfd594c0574b
md"#### Named tuples"

# ╔═╡ a57d4a77-7dd9-4a16-a853-1b15cc7fece2
md"###### Ejemplo"

# ╔═╡ d460cdee-f779-4463-8be2-6e9bfa20eb18
let
	p = (x = 1.1, y = 2.4) #NamedTuple{(:x, :y, Tuple{Float64, Float64}}

	K = keys(p)
	display(K)
	V = values(p)
	display(V)
end

# ╔═╡ 90823769-e873-4f42-9c8c-de580b897f6b
md"""
### Diccionarios

Un diccionario es una colección de pares de key-values
"""

# ╔═╡ ff5c343a-79e4-41da-bdfb-5407f174de3b
md"#### Creando un diccionario"

# ╔═╡ f6b20f93-b243-49b0-95ef-37c1f18cda4e
md"###### Ejemplo"

# ╔═╡ 4bbb9d60-5c3e-4eea-a801-6ae7f49d735e
let
	D1 = Dict("a" => 1, "b"=> 2, 1 => "a")
	display(D1)

	D2 = Dict([("a",1),("b",2),(1,"a")])
	display(D2)
end

# ╔═╡ c07d7369-34e5-427c-b0cd-191f14491b96
md"#### Accediendo a los elementos"

# ╔═╡ a694b6a4-291e-48ab-bc16-ce8473e328f7
md"###### Ejemplo"

# ╔═╡ 84e4a044-debd-4ab1-9a0a-f5ae6367b6f3
let
	D = Dict([("a", 1), ("b", 2), (1,"a")])
	#D["a"]
	D[1]
end

# ╔═╡ 9143f3d4-7ea1-4109-a6c3-090f5795ca3d
let
	D = Dict([("a", 1), ("b", 2), (1,"a")])
	for e in D
		println(e)
	end
end

# ╔═╡ 02250233-f7a5-4b42-94f5-098130665f76
let
	D = Dict([("a", 1), ("b", 2), (1,"a")])
	for (k,v) in D
		println(k,"=>",v)
	end
end

# ╔═╡ d4ed35fa-f148-4e96-9ce8-11f4dcdb8847
md"#### Modificando un diccionario"

# ╔═╡ 0249e76f-cb4c-4d7c-abd6-814eea00901c
md"###### Ejemplo"

# ╔═╡ ab0b151a-f727-4df3-8051-3eca6f9ae283
let
	D = Dict([("a", 1), ("b", 2), (1,"a")])
	D["c"] = 3
	display(D)
	D["c"] = "Hola"
	display(D)
	D = delete!(D,"c")
	display(D)
end

# ╔═╡ 28c02c9f-a88f-4859-8e47-223e91b327bf
md"""
### Estructuras
Las estructuras  definen un nuevo tipo coompuesto, están basadas en campos dados y tienen diferentes tipos. Una vez inicializadas no hay forma de modificarlas.
""" 

# ╔═╡ ed673bd2-404a-40d3-a533-34702aac012c
md"""
###### Ejemplo
```julia
struct Location
name::String
lat::Float64
lon::Float64
end
```
"""

# ╔═╡ 270dd99f-0b99-4254-95fc-d0237fb1a6e8
let
	struct Location
	name::String	
	lat::Float64
	lon::Float64
	end

	loc1 = Location("Los Angeles", 34.0522,-118.2437)

	display(loc1.name)
	display(loc1.lat)
	
end

# ╔═╡ 28bc0c9b-26ff-4e3d-9396-5bf88f1fd1cb
md"###### Ejemplo: Difinir vector con una estructura"

# ╔═╡ 6d23f861-2870-4cdf-8382-e4605aae4e4f
let
	sites = Location[]
	push!(sites,Location("Los Angeles", 34.0522,-118.2437))
	push!(sites,Location("Las Vegas", 36.1699,-115.1398))

	display(sites)
end

# ╔═╡ a03a2ecc-fead-4494-b243-db67fb4968a3
md"""
#### Estructuras mutables
Si queremos modificar las componentes de una estructura después de haber sido inicializada, debemos usar una estructrua mutable.
```julia
mutable struct mLocation
    name::String
    lat::Float32
    lon::Float32
end
```

"""

# ╔═╡ d5564a75-1b64-49fe-9335-f5701c90d4c9


# ╔═╡ 1af54509-5192-4646-a517-0b74f912e1a9
md"""
###### Ejemplo
"""

# ╔═╡ 22654a7b-0ed3-4917-aee8-a4e4c4e2f73f
let
	mutable struct mLocation
    	name::String
    	lat::Float32
    	lon::Float32
	end	
end

# ╔═╡ 689e6ccd-db5d-4c05-a1ef-d5a42b081a1a
let
loc1 = mLocation("Los Angeles", 34.0522,-118.2437)
loc1.name = "LA"
end

# ╔═╡ c32eeb2d-5a2c-4a9b-9f6d-406329c9020d


# ╔═╡ 977fc666-2229-4620-826e-5b92a3fa9dfe
md"""
## Gráficas
Existen varios paquetes para graficar
- Plots
- Plotly
- Pyplot
- Makie
- ...
"""

# ╔═╡ 2ee8573b-d524-4fd5-bfb8-7714efd57324
md"###### Ejemplos"

# ╔═╡ 9651bc7d-e0b7-4d1c-9933-bb1f7233ed4c
let
x = 0:0.05:1;
y = sin.(2π*x);
plot(x,y)
end

# ╔═╡ b7d41a4f-0767-4a9b-82fd-2e9c89ab9165
let
x = 0:0.05:1;
y = sin.(2π*x);
plot(x,y,
	seriestype =:path,
	linestyle = :dash
)
end

# ╔═╡ a96452c6-c876-476c-a867-33d0cea9fa8e
let
	x = 0:0.05:1;
y = sin.(2π*x);

plot(   x, y, 
        seriestype=:sticks,
        linestyle=:dash,
        lw = 3,
        seriescolor = :green,
        marker = :circle,
        markersize = 8,
        markercolor = :green,
        markerstrokecolor = :green,
    )
end

# ╔═╡ 64c1d447-f64e-47f2-874a-6d36c74e10af
let
x = 0:0.05:1;
y = sin.(2π*x);

plot(   x, y, 
        title="title (titlefontsize)",
        xlabel="xlabel (guidefontsize)",
        ylabel="ylabel (guidefontsize)",
        label="sin(2*pi*x) (legendfontsize)",
        titlefontsize=18,
        guidefontsize=18,
        tickfontsize=16,
        legendfontsize=12,
        legend=:topright
    )
end

# ╔═╡ 7ace162b-e7e4-41e4-8c5c-ff3fb8326654
let
	x = 0:0.05:1;
y = sin.(2π*x);
plot(   x, y, 
        line=(3,:green,:dash,:sticks),
        marker=(:circle,8,:green,:green),
        title="title",
        xlabel="xlabel",
        ylabel="ylabel",
        label=L"\sin(2\pi x)",
        legend=:outertopright,
        titlefontsize=18,
        guidefontsize=18,
        tickfontsize=16,
        legendfontsize=18,
        grid=false
    )
end

# ╔═╡ 76967e7d-1fcb-45a7-96de-0973c9865d4c
md"""
Algunas de las opciones para la posición de `legend` son:
- :right
- :left
- :top
- :bottom
- :topright
- :topleft
- :bottomright
- :bottomleft
- :outerright
- :outerleft
- :outertopright
- :outertopleft
- :outerbottomright
- :outerbottomleft
"""

# ╔═╡ 8c5db09c-f067-4c62-9c86-2a06cd89fc88
let
x = 0:0.05:1;
y = sin.(2π*x);
p1 = scatter(x, y, label="(x,y)")    
p2 = scatter(x, -y, label="(x,-y)", legend=:topleft)    
p3 = scatter(-x, y, label="(-x,y)", legend=:topleft)    
p4 = scatter(-x, -y, label="(-x,-y)")  

plot(p1,p2,p3,p4, layout = (2,2) )	
end

# ╔═╡ c10457ee-1b61-4103-bcde-9c5e0b415bcf
let
x = 0:0.01:1;
y = @. sin(8*2π*x)/exp(10*x);

plot(
    x,y, 
    lw=3,
    label=false,
    ylimits=(-1,1),
    xlimits=(0,0.5)
    )
end

# ╔═╡ 8ab71073-4511-4bc6-a373-17e6309435b2
let
	f(x) = exp(x)
p1 = plot(f, 1, 10, lw=3,
            title="Regular Plot of exp(x)",
            label=false)

p2 = plot(f, 1, 10, 
            yscale=:log10, lw=3,
            title="Semilog-y plot of exp(x)",
            label=false)

plot(p1,p2,layout=(2,1))
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Images = "916415d5-f1e6-5110-898d-aaa5f9f070e0"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
Images = "~0.25.2"
LaTeXStrings = "~1.3.0"
Plots = "~1.38.2"
PlutoUI = "~0.7.49"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.5"
manifest_format = "2.0"
project_hash = "9497016a148d33066ade0c80b77789feb15ebac9"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "69f7020bd72f069c219b5e8c236c1fa90d2cb409"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.2.1"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "195c5505521008abea5aee4f96930717958eac6f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.4.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "1dd4d9f5beebac0c03446918741b1a03dc5e5788"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.6"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.CatIndices]]
deps = ["CustomUnitRanges", "OffsetArrays"]
git-tree-sha1 = "a0f80a09780eed9b1d106a1bf62041c2efc995bc"
uuid = "aafaddc9-749c-510e-ac4f-586e18779b91"
version = "0.2.2"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c6d890a52d2c4d55d326439580c3b8d0875a77d9"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.7"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "38f7a08f19d8810338d4f5085211c7dfa5d5bdd8"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.4"

[[deps.Clustering]]
deps = ["Distances", "LinearAlgebra", "NearestNeighbors", "Printf", "Random", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "64df3da1d2a26f4de23871cd1b6482bb68092bd5"
uuid = "aaaa29a8-35af-508c-8bc3-b662a17a0fe5"
version = "0.14.3"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random", "SnoopPrecompile"]
git-tree-sha1 = "aa3edc8f8dea6cbfa176ee12f7c2fc82f0608ed3"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.20.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "00a2cccc7f098ff3b66806862d275ca3db9e6e5a"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.5.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.1+0"

[[deps.ComputationalResources]]
git-tree-sha1 = "52cb3ec90e8a8bea0e62e275ba577ad0f74821f7"
uuid = "ed09eef8-17a6-5b46-8889-db040fac31e3"
version = "0.3.2"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.CoordinateTransformations]]
deps = ["LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "681ea870b918e7cff7111da58791d7f718067a19"
uuid = "150eb455-5306-5404-9cee-2592286d6298"
version = "0.6.2"

[[deps.CustomUnitRanges]]
git-tree-sha1 = "1a3f97f907e6dd8983b744d2642651bb162a3f7a"
uuid = "dc8bdbbb-1ca9-579f-8c36-e416f6a65cce"
version = "1.0.2"

[[deps.DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FFTViews]]
deps = ["CustomUnitRanges", "FFTW"]
git-tree-sha1 = "cbdf14d1e8c7c8aacbe8b19862e0179fd08321c2"
uuid = "4f61f5a4-77b1-5117-aa51-3ab5ef4ef0cd"
version = "0.3.2"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "90630efff0894f8142308e334473eba54c433549"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.5.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "7be5f99f7d15578798f338f5433b6c432ea8037b"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "387d2b8b3ca57b791633f0993b31d8cb43ea3292"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.71.3"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "5982b5e20f97bff955e9a2343a14da96a746cd8c"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.71.3+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "78e2c69783c9753a91cdae88a8d432be85a2ab5e"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "d3b3624125c1474292d0d8ed0f65554ac37ddb23"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+2"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "d61890399bc535850c4bf08e4e0d3a7ad0f21cbd"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "ba2d094a88b6b287bd25cfa86f301e7693ffae2f"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.7.4"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "Dates", "IniFile", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "eb5aa5e3b500e191763d35198f859e4b40fff4a6"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.7.3"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "c54b581a83008dc7f292e205f4c409ab5caa0f04"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.10"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "b51bb8cae22c66d0f6357e3bcb6363145ef20835"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.5"

[[deps.ImageContrastAdjustment]]
deps = ["ImageCore", "ImageTransformations", "Parameters"]
git-tree-sha1 = "0d75cafa80cf22026cea21a8e6cf965295003edc"
uuid = "f332f351-ec65-5f6a-b3d1-319c6670881a"
version = "0.3.10"

[[deps.ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "acf614720ef026d38400b3817614c45882d75500"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.9.4"

[[deps.ImageDistances]]
deps = ["Distances", "ImageCore", "ImageMorphology", "LinearAlgebra", "Statistics"]
git-tree-sha1 = "b1798a4a6b9aafb530f8f0c4a7b2eb5501e2f2a3"
uuid = "51556ac3-7006-55f5-8cb3-34580c88182d"
version = "0.2.16"

[[deps.ImageFiltering]]
deps = ["CatIndices", "ComputationalResources", "DataStructures", "FFTViews", "FFTW", "ImageBase", "ImageCore", "LinearAlgebra", "OffsetArrays", "Reexport", "SnoopPrecompile", "SparseArrays", "StaticArrays", "Statistics", "TiledIteration"]
git-tree-sha1 = "f265e53558fbbf23e0d54e4fab7106c0f2a9e576"
uuid = "6a3955dd-da59-5b1f-98d4-e7296123deb5"
version = "0.7.3"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "342f789fd041a55166764c351da1710db97ce0e0"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.6"

[[deps.ImageMagick]]
deps = ["FileIO", "ImageCore", "ImageMagick_jll", "InteractiveUtils", "Libdl", "Pkg", "Random"]
git-tree-sha1 = "5bc1cb62e0c5f1005868358db0692c994c3a13c6"
uuid = "6218d12a-5da1-5696-b52f-db25d2ecc6d1"
version = "1.2.1"

[[deps.ImageMagick_jll]]
deps = ["Artifacts", "Ghostscript_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "124626988534986113cfd876e3093e4a03890f58"
uuid = "c73af94c-d91f-53ed-93a7-00f77d67a9d7"
version = "6.9.12+3"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "36cbaebed194b292590cba2593da27b34763804a"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.8"

[[deps.ImageMorphology]]
deps = ["ImageCore", "LinearAlgebra", "Requires", "TiledIteration"]
git-tree-sha1 = "e7c68ab3df4a75511ba33fc5d8d9098007b579a8"
uuid = "787d08f9-d448-5407-9aad-5290dd7ab264"
version = "0.3.2"

[[deps.ImageQualityIndexes]]
deps = ["ImageContrastAdjustment", "ImageCore", "ImageDistances", "ImageFiltering", "LazyModules", "OffsetArrays", "SnoopPrecompile", "Statistics"]
git-tree-sha1 = "5985d467623f106523ed8351f255642b5141e7be"
uuid = "2996bd0c-7a13-11e9-2da2-2f5ce47296a9"
version = "0.3.4"

[[deps.ImageSegmentation]]
deps = ["Clustering", "DataStructures", "Distances", "Graphs", "ImageCore", "ImageFiltering", "ImageMorphology", "LinearAlgebra", "MetaGraphs", "RegionTrees", "SimpleWeightedGraphs", "StaticArrays", "Statistics"]
git-tree-sha1 = "36832067ea220818d105d718527d6ed02385bf22"
uuid = "80713f31-8817-5129-9cf8-209ff8fb23e1"
version = "1.7.0"

[[deps.ImageShow]]
deps = ["Base64", "FileIO", "ImageBase", "ImageCore", "OffsetArrays", "StackViews"]
git-tree-sha1 = "b563cf9ae75a635592fc73d3eb78b86220e55bd8"
uuid = "4e3cecfd-b093-5904-9786-8bbb286a6a31"
version = "0.3.6"

[[deps.ImageTransformations]]
deps = ["AxisAlgorithms", "ColorVectorSpace", "CoordinateTransformations", "ImageBase", "ImageCore", "Interpolations", "OffsetArrays", "Rotations", "StaticArrays"]
git-tree-sha1 = "8717482f4a2108c9358e5c3ca903d3a6113badc9"
uuid = "02fcd773-0e25-5acc-982a-7f6622650795"
version = "0.9.5"

[[deps.Images]]
deps = ["Base64", "FileIO", "Graphics", "ImageAxes", "ImageBase", "ImageContrastAdjustment", "ImageCore", "ImageDistances", "ImageFiltering", "ImageIO", "ImageMagick", "ImageMetadata", "ImageMorphology", "ImageQualityIndexes", "ImageSegmentation", "ImageShow", "ImageTransformations", "IndirectArrays", "IntegralArrays", "Random", "Reexport", "SparseArrays", "StaticArrays", "Statistics", "StatsBase", "TiledIteration"]
git-tree-sha1 = "03d1301b7ec885b266c0f816f338368c6c0b81bd"
uuid = "916415d5-f1e6-5110-898d-aaa5f9f070e0"
version = "0.25.2"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "87f7662e03a649cffa2e05bf19c303e168732d3e"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.2+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "5cd07aab533df5170988219191dfad0519391428"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.3"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.IntegralArrays]]
deps = ["ColorTypes", "FixedPointNumbers", "IntervalSets"]
git-tree-sha1 = "be8e690c3973443bec584db3346ddc904d4884eb"
uuid = "1d092043-8f09-5a30-832f-7509e371ab51"
version = "0.1.5"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "721ec2cf720536ad005cb38f50dbba7b02419a15"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.14.7"

[[deps.IntervalSets]]
deps = ["Dates", "Random", "Statistics"]
git-tree-sha1 = "16c0cc91853084cb5f58a78bd209513900206ce6"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.4"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.JLD2]]
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "Pkg", "Printf", "Reexport", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "ec8a9c9f0ecb1c687e34c1fda2699de4d054672a"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.29"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "f377670cda23b6b7c1c0b3893e37451c5c1a2185"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.5"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "a77b273f1ddec645d1b7c4fd5fb98c8f90ad10a5"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.1"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "2422f47b34d4b127720a18f86fa7b1aa2e141f29"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.18"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "946607f84feb96220f480e0422d3484c49c00239"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.19"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "cedb76b37bc5a6c702ade66be44f831fa23c681e"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "2ce8695e1e699b68702c03402672a69f54b8aca9"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2022.2.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.MappedArrays]]
git-tree-sha1 = "e8b359ef06ec72e8c030463fe02efe5527ee5142"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.1"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.MetaGraphs]]
deps = ["Graphs", "JLD2", "Random"]
git-tree-sha1 = "1130dbe1d5276cb656f6e1094ce97466ed700e5a"
uuid = "626554b9-1ddb-594c-aa3c-2596fe9399a5"
version = "0.7.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "a7c3d1da1189a1c2fe843a3bfa04d18d20eb3211"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.1"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "2c3726ceb3388917602169bed973dbc97f1b51a8"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.13"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "5ae7ca23e13855b3aba94550f26146c01d259267"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.0"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "f71d8950b724e9ff6110fc948dff5a329f901d64"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.8"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "923319661e9a22712f24596ce81c54fc0366f304"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.1.1+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "6503b77492fd7fcb9379bf73cd31035670e3c509"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.3.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6e9dba33f9f2c44e08a020b0caf6903be540004"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.19+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.40.0+0"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "f809158b27eba0c18c269cf2a2be6ed751d3e81d"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.3.17"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "03a7a85b76381a3d04c7a1656039197e70eda03d"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.11"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "8175fc2b118a3755113c8e68084dc1a9e63c61ee"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.3"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f6cf8e7944e50901594838951729a1861e668cb8"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.2"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "SnoopPrecompile", "Statistics"]
git-tree-sha1 = "5b7690dd212e026bbab1860016a6601cb077ab66"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.2"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SnoopPrecompile", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "a99bbd3664bb12a775cda2eba7f3b2facf3dad94"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.38.2"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "eadad7b14cf046de6eb41f13c9275e5aa2711ab6"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.49"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "d7a7aef8f8f2d537104f170139553b14dfe39fe9"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.2"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "0c03844e2231e12fda4d0086fd7cbe4098ee8dc5"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+2"

[[deps.Quaternions]]
deps = ["LinearAlgebra", "Random", "RealDot"]
git-tree-sha1 = "da095158bdc8eaccb7890f9884048555ab771019"
uuid = "94ee1d12-ae83-5a48-8b1c-48b8ff168ae0"
version = "0.7.4"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "dc84268fe0e3335a62e315a3a7cf2afa7178a734"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.3"

[[deps.RealDot]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9f0a1b71baaf7650f4fa8a1d168c7fb6ee41f0c9"
uuid = "c1ae055f-0cd5-4b69-90a6-9a35b1a98df9"
version = "0.1.0"

[[deps.RecipesBase]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "261dddd3b862bd2c940cf6ca4d1c8fe593e457c8"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.3"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase", "SnoopPrecompile"]
git-tree-sha1 = "e974477be88cb5e3040009f3767611bc6357846f"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.11"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RegionTrees]]
deps = ["IterTools", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "4618ed0da7a251c7f92e869ae1a19c74a7d2a7f9"
uuid = "dee08c22-ab7f-5625-9660-a9af2021b33f"
version = "0.3.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rotations]]
deps = ["LinearAlgebra", "Quaternions", "Random", "StaticArrays", "Statistics"]
git-tree-sha1 = "9480500060044fd25a1c341da53f34df7443c2f2"
uuid = "6038ab10-8711-5258-84ad-4b1120ba62dc"
version = "1.3.4"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "f94f779c94e58bf9ea243e77a37e16d9de9126bd"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SimpleWeightedGraphs]]
deps = ["Graphs", "LinearAlgebra", "Markdown", "SparseArrays", "Test"]
git-tree-sha1 = "a6f404cc44d3d3b28c793ec0eb59af709d827e4e"
uuid = "47aef6b3-ad0c-573a-a1e2-d07658019622"
version = "1.2.1"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "8fb59825be681d451c246a795117f317ecbcaa28"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.2"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "6954a456979f23d05085727adb17c4551c19ecd1"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.12"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "UUIDs"]
git-tree-sha1 = "7e6b0e3e571be0b4dd4d2a9a3a83b65c04351ccc"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.6.3"

[[deps.TiledIteration]]
deps = ["OffsetArrays"]
git-tree-sha1 = "5683455224ba92ef59db72d10690690f4a8dc297"
uuid = "06e1c1a7-607b-532d-9fad-de7d9aa2abac"
version = "0.3.1"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "94f38103c984f89cf77c402f2a68dbd870f8165f"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.11"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.URIs]]
git-tree-sha1 = "ac00576f90d8a259f2c9d823e91d1de3fd44d348"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "ed8d92d9774b077c53e1da50fd81a36af3744c1c"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "93c41695bc1c08c46c5899f4fe06d6ead504bb73"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.3+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "868e669ccb12ba16eaf50cb2957ee2ff61261c56"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.29.0+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "libpng_jll"]
git-tree-sha1 = "d4f63314c8aa1e48cd22aa0c17ed76cd1ae48c3c"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.3+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9ebfc140cc56e8c2156a15ceac2f0302e327ac0a"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+0"
"""

# ╔═╡ Cell order:
# ╠═a22dba1a-9881-11ed-3f52-559778073a9c
# ╟─df6f74e1-42eb-4e28-9b7f-1d8bb8bb5b0f
# ╟─331513e1-d6c0-4f4f-819b-4e32319bce80
# ╟─b11c332f-2375-45f7-a4ab-72ed58a5c3f8
# ╟─af48b15a-f04c-42ab-a1fe-8ff965e0a6d4
# ╟─e5a72dd2-2d29-42dc-b125-93e7072c2697
# ╟─802e3512-814d-48b7-b5b6-54ecbdcf63bb
# ╟─95c959ab-9900-4290-98ac-57af9bec14a9
# ╟─754d000d-1f35-4a29-ad1e-0ab6316e33a1
# ╟─0966daf4-8e14-4569-b68a-54f24dbf04e5
# ╟─706b8cc1-5b96-4bda-ae4f-fd32d6842187
# ╟─b612721a-1b06-49a8-86b7-8a3f09480597
# ╟─0edc622e-cdd4-4595-a43b-eed80ab301a8
# ╟─5f1e7909-583b-42d1-a9bb-9e63ae7ceb99
# ╟─ed620ff7-9458-4d8e-9d51-62fe0a97c623
# ╟─6d2c74a1-3fd1-4310-b6a1-c2eac4243257
# ╟─00fe7ab2-b468-4489-aade-9935d377d236
# ╟─83fa2bd9-b168-46c2-87d2-2d4fa2844e47
# ╟─0adba4a5-c418-43b5-aa08-6b29aaf025c1
# ╟─06f24748-40fb-480f-a202-4de665de2ec8
# ╟─65bc4440-574c-409b-95ce-970584670d42
# ╟─02d484f5-0090-42a8-8c89-7f76b19d4ac5
# ╟─d8ee737b-3e56-4f71-9143-cbef1c4ec618
# ╟─9741bd55-5d21-4a1e-9f0f-2819de092ae2
# ╟─b3d8fcab-9413-43cd-a032-80f837dd6896
# ╟─82c606ec-06ea-4ab8-89e3-a84e314bc5bd
# ╟─e929b471-563c-49db-8c6d-01b5924e4e03
# ╟─9be54d29-b85e-4337-bf0a-bb2cc1e13e02
# ╟─82dfe60b-2779-4ec0-96eb-023b9a1d84d1
# ╟─63573df9-2133-4e10-89b6-b33a7fe5d3de
# ╟─f11eaf6a-e3d0-4ed5-b64e-1b0a27ffe398
# ╟─e3e66c98-c462-4cb9-86f0-d462a005977b
# ╟─7f01e0e9-9885-4d85-802a-457980507d87
# ╟─f8bca9af-a155-496d-a1c7-2eb0b0822037
# ╟─ed3951f4-df95-468f-94b9-c694d13c1605
# ╟─ffd77c4b-99c6-4591-95eb-15ca69ed5361
# ╟─eaadf6cf-0cf5-43d9-943d-3975fc8a356d
# ╟─182e369d-b31b-4f7b-b455-55cdfc0fb4b1
# ╟─4a094840-92f5-4447-aa9d-42690f40dcb2
# ╟─ee262cee-ab01-49cb-bc70-2c73af12f1cd
# ╟─ab3980f9-a467-4212-8c18-dd1b2442f07c
# ╟─c944d642-b67a-41fa-a193-c156fdc02479
# ╟─87f969e4-e208-4dcf-b548-a1684345969e
# ╟─02aa2202-5ec1-404e-a7bf-33b974c47588
# ╟─f39e45d1-33f3-4a5a-b152-5f39be2d6f89
# ╟─d8928380-e9d8-4004-a8d2-fae26cbe6a6e
# ╟─05ee5b6d-b6f0-43eb-a794-ad5131d6b1f4
# ╟─c12fcbbd-46d3-413b-aff2-7631ae108cc4
# ╟─d5b049b8-15ac-46a0-9145-bd3ea24c8aac
# ╟─3c3b5816-928c-4c00-959d-6c33da5ec111
# ╟─0360935c-09dc-4775-94bb-bc54dd631138
# ╟─07f1cba5-eae0-4b5f-8c62-3557a71accc1
# ╟─7d440a4f-2a47-453f-a47f-15a67964cc36
# ╟─a9a5b3f6-2019-4989-b04c-f55d2437f530
# ╠═920bfd76-525f-4a5f-b4f5-2a330e620c8c
# ╟─3fcee3a8-14ad-4dba-a6b6-b729ae0639e4
# ╟─0d8b4d11-664e-4cd1-8d24-704870f65a1c
# ╟─040e1ac8-b912-4e2b-a5de-8053133a0882
# ╟─91e114d7-5347-462e-a06c-dd05fcb106c7
# ╟─a9745a65-c58b-49ab-b47d-a1df7466cad1
# ╟─2a79ca9b-dc5a-4198-990e-69967d0e3c82
# ╟─3b06956a-622f-424c-995c-679425c354c0
# ╟─71362aea-97eb-4fbf-b5f7-cb3644a70045
# ╟─662e3bb5-32ac-4f79-ba9a-51465b420a9e
# ╟─fa4a472f-0972-46ba-b31c-51bd9ba23ff2
# ╟─d9cdaec4-91ac-4105-93f9-d6c393a3fa73
# ╟─e9d1af7e-c8d4-4231-a600-373b56d7cefe
# ╟─18479558-da15-4343-a51c-83aa20e92eed
# ╟─8eab62ad-4cc3-4b76-b02e-2f5b942be980
# ╟─f70276aa-3e9c-42f8-9021-46fccdf3802c
# ╟─b8883753-a281-4c56-932e-b39987a67f97
# ╠═1abf3f43-9f3d-407d-ab26-3aef2efe5d58
# ╟─52cf0faf-4bff-4f55-88a1-e684810cdf03
# ╠═34613b4d-34fa-4687-89f0-79d8a150486c
# ╟─24d42c25-3457-43c3-8ff6-d74a7719113a
# ╟─678268a0-9dc2-4f5f-a3fb-6ee906401aef
# ╠═1f2a56c3-f5de-453c-b965-9c066fd67aaa
# ╟─b48aac04-5212-4f11-9884-b45564157174
# ╠═04950b12-3b15-4ac5-96f6-18963a3814a0
# ╟─8ec133a1-76ac-4bca-924e-6880dfc6633c
# ╠═021e59a9-bd36-46e1-aa79-67e3a1953dfe
# ╟─ff67f2aa-f208-48a9-8422-0d795145af62
# ╠═bde8900d-6907-4b2c-b0f6-44428d83a684
# ╟─11d479a7-71a2-4879-bc16-dc6053cb8c76
# ╟─fdfcd417-1d9c-4fd4-9df8-d50a59bddc1b
# ╟─5bd0a998-bcbd-4543-a899-480d0be845a7
# ╠═e51bd7c6-8649-4140-be0c-ae66902e834a
# ╠═6faf62b9-003d-487e-98fe-23f5940c7562
# ╟─7abe3e0f-6b6d-47a4-a86b-a7b666715ac2
# ╠═ac884bda-f2a2-4666-bbd5-76339087dde3
# ╟─34ec3b4d-2d61-4666-ab76-9d71c274d582
# ╠═980a69af-e894-4d19-b5bf-9fac8f2900f5
# ╟─27b65217-9841-4d6e-b60e-345b0c5e1245
# ╠═25585d86-1909-47f7-92b8-297ef9b8c10b
# ╟─b2b54e1f-89b9-4f0b-9d73-3d14dd6a7a2a
# ╟─a65cb56f-0e37-47d3-939b-6b72aa23f230
# ╟─c9b4affa-1e71-4b44-a92e-d914013e8747
# ╠═f9119e93-da01-488b-817b-f8656e3b72b8
# ╟─0cf416ca-ecc0-4193-9c03-18bdfb50d737
# ╠═f4b3fec0-6786-415a-b709-9f311c50a6b2
# ╠═aca8a8ff-2708-4ccc-b6c7-027697bc5b7f
# ╟─94842940-9d38-4a36-bf26-382fa5ce15ef
# ╠═0e8ab455-8bdd-4688-b483-3e8cd3cb0b8f
# ╟─a22a5694-672b-489d-a272-b4aafa7c2d9d
# ╟─89813b0d-1884-4888-b540-fa938d8ac0f4
# ╠═55d97f3b-873f-410d-b245-a4222072d6b6
# ╟─bc1e6fbc-4a1a-47cf-ae79-840fad31c1fb
# ╟─0628ca0c-1b8e-4f01-aa79-e7c99bc8f5f7
# ╟─65449b04-9873-4c1d-aa8a-3a3133823424
# ╠═26736f5a-64b8-4f39-a8cd-71218a922726
# ╠═fbf289d9-e3ca-43c1-96ad-e33f35addbf5
# ╟─3a129775-7338-4c2b-b1de-92cf3ea5897e
# ╟─c1c33a9b-43be-4257-b44f-ada3bbfc3f50
# ╠═47be421e-24c0-453c-a8b6-eb9a118a26ec
# ╟─552daa2c-4ff9-4d03-b697-bfd594c0574b
# ╟─a57d4a77-7dd9-4a16-a853-1b15cc7fece2
# ╠═d460cdee-f779-4463-8be2-6e9bfa20eb18
# ╟─90823769-e873-4f42-9c8c-de580b897f6b
# ╟─ff5c343a-79e4-41da-bdfb-5407f174de3b
# ╟─f6b20f93-b243-49b0-95ef-37c1f18cda4e
# ╠═4bbb9d60-5c3e-4eea-a801-6ae7f49d735e
# ╟─c07d7369-34e5-427c-b0cd-191f14491b96
# ╟─a694b6a4-291e-48ab-bc16-ce8473e328f7
# ╠═84e4a044-debd-4ab1-9a0a-f5ae6367b6f3
# ╠═9143f3d4-7ea1-4109-a6c3-090f5795ca3d
# ╠═02250233-f7a5-4b42-94f5-098130665f76
# ╟─d4ed35fa-f148-4e96-9ce8-11f4dcdb8847
# ╟─0249e76f-cb4c-4d7c-abd6-814eea00901c
# ╠═ab0b151a-f727-4df3-8051-3eca6f9ae283
# ╟─28c02c9f-a88f-4859-8e47-223e91b327bf
# ╟─ed673bd2-404a-40d3-a533-34702aac012c
# ╠═270dd99f-0b99-4254-95fc-d0237fb1a6e8
# ╟─28bc0c9b-26ff-4e3d-9396-5bf88f1fd1cb
# ╠═6d23f861-2870-4cdf-8382-e4605aae4e4f
# ╟─a03a2ecc-fead-4494-b243-db67fb4968a3
# ╟─d5564a75-1b64-49fe-9335-f5701c90d4c9
# ╟─1af54509-5192-4646-a517-0b74f912e1a9
# ╠═22654a7b-0ed3-4917-aee8-a4e4c4e2f73f
# ╠═689e6ccd-db5d-4c05-a1ef-d5a42b081a1a
# ╠═c32eeb2d-5a2c-4a9b-9f6d-406329c9020d
# ╟─977fc666-2229-4620-826e-5b92a3fa9dfe
# ╠═e029e588-7e17-4888-bb14-a7aa94a95f82
# ╟─2ee8573b-d524-4fd5-bfb8-7714efd57324
# ╠═9651bc7d-e0b7-4d1c-9933-bb1f7233ed4c
# ╠═b7d41a4f-0767-4a9b-82fd-2e9c89ab9165
# ╠═a96452c6-c876-476c-a867-33d0cea9fa8e
# ╠═64c1d447-f64e-47f2-874a-6d36c74e10af
# ╠═066a6395-81cb-4323-85da-31deccc9f7e0
# ╠═7ace162b-e7e4-41e4-8c5c-ff3fb8326654
# ╟─76967e7d-1fcb-45a7-96de-0973c9865d4c
# ╠═8c5db09c-f067-4c62-9c86-2a06cd89fc88
# ╠═c10457ee-1b61-4103-bcde-9c5e0b415bcf
# ╠═8ab71073-4511-4bc6-a373-17e6309435b2
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002