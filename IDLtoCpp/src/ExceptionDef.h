/**
 *  @file ExceptionDef.h
 *  @date Created on: Apr 12, 2013
 *  @author Sabatino Severino (sabatino.severino@remocean.com)
 *  @bug No known bugs.
 */

#ifndef EXCEPTIONDEF_H_
#define EXCEPTIONDEF_H_

//
// Eccezioni
//

/** @def EXCP_UNHANDLED Eccezione non ancora gestita */
#define EXCP_UNHANDLED 0

/** @def EXCP_UNHANDLED Result dimensions must be integer factor of original dimensions */
#define EXCP_REBIN_DIM_NOT_VALID 1

#define EXCP_SUM_BOUND_NOT_VALID 2
#define EXCP_MINUS_BOUND_NOT_VALID 3
#define EXCP_MULT_BOUND_NOT_VALID 4

/**
 * brief Converte l'enum relativo all'eccezione in una corrispettiva String
 * @param [in] e Tipo di eccezione da convertire
 * @return Eccezione e convertita in String
 */
static std::string ExceptionMsg(int e) {

	std::string msg;
	switch (e) {
	case (EXCP_UNHANDLED):
		msg = "Eccezione non gestita";
		break;
	case (EXCP_REBIN_DIM_NOT_VALID):
		msg = "Result dimensions must be integer factor of original dimensions";
		break;
	case (EXCP_SUM_BOUND_NOT_VALID):
		msg = "Sum operator: Dimensions not valid";
		break;
	case (EXCP_MINUS_BOUND_NOT_VALID):
		msg = "Minus operator: Dimensions not valid";
		break;
	case (EXCP_MULT_BOUND_NOT_VALID):
		msg = "Mult operator: Dimensions not valid";
		break;

	default:
		msg = "EXCP_UNKNOW";
	}
	return msg;
}

#endif /* EXCEPTIONDEF_H_ */
